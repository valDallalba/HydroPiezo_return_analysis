
import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import folium
import io
from pyproj import Transformer
import requests

####
#Function to analyse the return period of piezometric levels using Gumbel distribution
####

def extract_annual_maxima(data):
    # Extraire les maxima annuels (valeurs maximales par année)
    extrem_event = data.loc[data.groupby(data.index.year)['nvx_piezo'].idxmax()]
    
    # Ajouter le mois du maximum pour information (facultatif)
    extrem_event['Month'] = extrem_event.index.month
    
    return extrem_event

def calculate_ranks(extrem_event):
    # Trier les données par ordre croissant
    donnees_triees = np.sort(extrem_event['nvx_piezo'].values)
    
    # Calculer les rangs de chaque donnée triée
    rangs = np.arange(1, len(donnees_triees) + 1)
    
    return donnees_triees, rangs

def calculate_hazen_frequenc_empirique(rangs, nb_data):
    # Calcul de la probabilité de Hazen pour chaque rang
    prob_hazen = (rangs - 0.5) / nb_data
    
    return prob_hazen

def calculate_reduced_variable_empirique(prob_hazen):
    # Calcul de la variable réduite u (Gumbel)
    u_gumbel = -np.log(-np.log(prob_hazen))
    
    return u_gumbel

def calculate_polyfit(rangs, donnees_triees, verbose=False):
    # Calcul du polyfit (ajustement linéaire) sur les données
    coefficients = np.polyfit(rangs, donnees_triees, 1)  # Ajustement linéaire (degré 1)
    slope, intercept = coefficients  # Extraire la pente (slope) et l'ordonnée à l'origine (intercept)
    
    # Affichage des coefficients
    if verbose:
        print(f"Coefficients du polyfit : Pente = {slope}, Ordonnée à l'origine = {intercept}")
    return coefficients

def value_from_return_period(intercep, slope, T):
    # Calculer la probabilité associée au temps de retour T
    # Le calcul differ pour les temps de retour court, la première formule est corrigé de RS
    # La seconde est la formule classique
    F_x = 1 + (1/(1-T))
    #F_x = 1 - (1 / T)  

    # Calculer la variable réduite u à partir de F(x)
    u = -np.log(-np.log(F_x))
    
    # Calculer la valeur x associée à u à partir des coefficients intercep (a) et slope (b)
    x = intercep + slope * u
    return x

def return_period_from_value(intercep, slope, x):
    # Calculer la variable réduite u à partir de la valeur x
    u = (x - intercep) / slope
    
    # Calculer la probabilité F(x) à partir de u
    F_x = np.exp(-np.exp(-u))
    
    # Calculer le temps de retour à partir de F(x)
    # Le calcul differ pour les temps de retour court, la première formule est corrigé de RS
    # La seconde est la formule classique
    T = 1 - (1/(F_x - 1))
    #T = 1 / (1 - F_x)
    return T

def max_serie_analysis(data):
    '''
    Analyse des séries maximales annuelles de niveaux piézométriques
    Parameters: data serie piezo
    Returns: DataFrame des résultats et liste des résultats intermédiaires
    1. data: série temporelle des niveaux piézométriques    
    2. extrem_event: maxima annuels extraits
    3. u_gumbel: variable réduite u (Gumbel)
    4. donnees_triees: données triées des maxima annuels
    5. rangs: rangs des données triées
    6. prob_hazen: probabilité empirique de Hazen
    7. slope: pente du polyfit
    8. intercept: ordonnée à l'origine du polyfit
    9. temps_retour: tableau des temps de retour estimés
    10. niveaux_retour: niveaux de retour estimés pour les temps de retour
    11. days_above_Q95: nombre de jours par année avec niveau piézométrique > Q95
    12. max_head_diff_per_year: différence maximale de niveau par année au-dessus de Q95
    '''
    # Étape 1: Extraire les maxima annuels
    extrem_event = extract_annual_maxima(data)
    
    # Étape 2: Calculer les rangs
    donnees_triees, rangs = calculate_ranks(extrem_event)
    
    # Étape 3: Calculer la fréquence empirique de Hazen
    nb_data    = len(donnees_triees)
    prob_hazen = calculate_hazen_frequenc_empirique(rangs, nb_data)
    
    # Étape 4: Calculer la variable réduite u (Gumbel)
    u_gumbel = calculate_reduced_variable_empirique(prob_hazen)
    
    # Création du DataFrame pour afficher les résultats
    df = pd.DataFrame({
        "Valeur Observée": donnees_triees,
        "Rang": rangs,
        "Probabilité Hazen": prob_hazen,
        "Variable réduite u (Gumbel)": u_gumbel
    })
    
    slope, intercept = calculate_polyfit(u_gumbel, donnees_triees, verbose=False)
    
    # Estimate levels for different return periods (10, 25, 50, 100 years)
    temps_retour   = np.arange(2.01, 200, 1)
    niveaux_retour = value_from_return_period(intercept, slope, temps_retour)

    # Calculate Q95 (95th percentile of water levels)
    Q95 = data['nvx_piezo'].quantile(0.95)

    # Calculate the number of days per year with water level above Q95
    days_above_Q95 = data[data['nvx_piezo'] > Q95].groupby(data[data['nvx_piezo'] > Q95].index.year).size()

    # For each year where there is a number of days above Q95, calculate the maximum head difference for the year
    max_head_diff_per_year = {}
    for year, days in days_above_Q95.items():
        if days > 0:
            yearly_data   = data[data.index.year == year]
            max_head_diff = yearly_data['nvx_piezo'].max() - Q95
            max_head_diff_per_year[year] = max_head_diff

    # Retourner le DataFrame et les résultats intermédiaires
    return df, [data, extrem_event, u_gumbel, donnees_triees, slope, intercept, temps_retour, niveaux_retour, days_above_Q95, max_head_diff_per_year]


######
# Functions for map creation and pdf creation
######

def get_high_res_swissalti3d_elevation(x, y):
    """
    Fetches the high-resolution altitude from SwissALTI3D using CH1903+ (LV95) coordinates.
    """
    url     = "https://api3.geo.admin.ch/rest/services/height"
    params  = {"easting": x,  # CH1903+ Easting (X)
               "northing": y,  # CH1903+ Northing (Y)
               "sr": 2056,  # CH1903+ / LV95
               "model": "swissALTI3D"  # High-resolution DEM model
        }
    
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        data      = response.json()
        elevation = data.get('height', 'N/A')
        return elevation
    
    else:
        return f"Error: {response.status_code}, {response.text}"


def create_map(x, y, zoom_start=17):
    """
    Create a map centered at the given coordinates with a specified zoom level.

    Parameters:
    x (float): Latitude coordinate in CH1903+.
    y (float): Longitude coordinate in CH1903+.
    zoom_start (int): Initial zoom level for the map. Default is 15.

    Returns:
    folium.Map: A folium map object centered at the given coordinates.
    """
    # Convert CH1903+ coordinates to WGS84
    transformer = Transformer.from_crs("epsg:2056", "epsg:4326")
    lat, lon    = transformer.transform(x, y)
    
    # Create a map centered at the given coordinates using Swiss GeoAdmin map
    map_ = folium.Map(location=[lat, lon], zoom_start=zoom_start, tiles='https://wmts.geo.admin.ch/1.0.0/ch.swisstopo.pixelkarte-farbe/default/current/3857/{z}/{x}/{y}.jpeg', attr='Swiss GeoAdmin')
    
    # Add a marker at the given coordinates with elevation info
    folium.Marker([lat, lon], icon=folium.Icon(color='red', icon_size=(40, 40))).add_to(map_)
    
    return map_


def save_map_as_image(map_, file_path):
    """
    Save a folium map as an image file.

    Parameters:
    map_ (folium.Map): The folium map object to save.
    file_path (str): The path to save the image file.
    """
    img_data = map_._to_png(5)
    img_ = Image.open(io.BytesIO(img_data))
    img_.save(file_path)
    return img_


def plot_max_serie_analysis(data_plot, name, path_to_save=False, to_plot=True):
    data, extrem_event, u_gumbel, donnees_triees, slope, intercept, temps_retour, niveaux_retour, days_above_Q95, max_head_diff_per_year = data_plot
    
    # Plot the full time series
    fig = plt.figure(figsize=(16, 12))  # A4 size in inches (21cm x 29.7cm)
    ax1 = fig.add_subplot(311)
    ax1.plot(data.index, data['nvx_piezo'], label="Série Temporelle", color="k", alpha=0.6)
    ax1.scatter(extrem_event.index, extrem_event['nvx_piezo'], color="darkred", label="Max Annuel", zorder=3)
    ax1.set_title("Série temporelle piézométrique et maximums annuels piézo {}".format(name))
    ax1.set_xlabel("Date")
    ax1.set_ylabel("Niveau piézométrique")

    # Add dashed lines for return periods
    return_periods = np.array([10, 25, 50, 100])
    niveaux        = value_from_return_period(intercept, slope, return_periods)
    colors         = list(mcolors.LinearSegmentedColormap.from_list("", ["blue", "purple"])(np.linspace(0, 1, len(return_periods))))

    for rp, nvx, color in zip(return_periods, niveaux, colors):
        ax1.axhline(y=nvx, color=color, linestyle='--', linewidth=0.8, label=f'Temps de retour {rp} ans: {nvx:.2f} m')

    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    ax1.grid(True, zorder=-4)

    # 2. Ajustement Gumbel (x vs u)
    ax2 = fig.add_subplot(323)
    ax2.scatter(u_gumbel, donnees_triees, label="Données", color='darkred')
    ax2.plot(u_gumbel, intercept + slope * u_gumbel, color="k", label=f"Régression: y = {slope:.2f}*u + {intercept:.2f}", linestyle='--')
    ax2.set_xlabel("Variable réduite u (Gumbel)")
    ax2.set_ylabel("Valeur Observée (m)")
    ax2.set_title("Ajustement des niveaux avec la loi de Gumbel")
    ax2.legend(loc='lower right')
    ax2.grid(True, zorder=-4)

    # 1. Distribution des valeurs par mois
    ax3 = fig.add_subplot(324)
    ax3.hist(extrem_event['Month'], bins=np.arange(1, 14) - 0.5, edgecolor='black', color='royalblue', alpha=0.95, width=0.8, zorder=4)
    ax3.set_xlabel("Mois")
    ax3.set_ylabel("Nombre d'événements")
    ax3.set_title("Distribution des maximums annuels par mois")
    ax3.set_xticks(np.arange(1, 13))
    ax3.set_xticklabels(['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Juin', 'Juil', 'Août', 'Sep', 'Oct', 'Nov', 'Déc'])
    ax3.grid(True, zorder=-4)

    # 3. Temps de retour vs Niveaux de retour
    ax4 = fig.add_subplot(325)
    ax4.plot(temps_retour, niveaux_retour, label="Courbe de retour", color='k')
    important_events = extrem_event.sort_values('nvx_piezo')[-3:]
    ax4.scatter(return_period_from_value(intercept, slope, important_events['nvx_piezo'].values), important_events['nvx_piezo'], color='darkred', zorder=5, label="3 évènements les plus importants")
    
    # Annotate the scatter plot with the year of the events
    for i, event in important_events.iterrows():
        ax4.annotate(i.year, (return_period_from_value(intercept, slope, event['nvx_piezo']), event['nvx_piezo']), 
                     textcoords="offset points", xytext=(0,10), ha='center', 
                     bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    # 4 . Add horizontal line for TN (Terrain Naturel)
    ax4.set_xlabel("Temps de retour [années]")
    ax4.set_ylabel("Niveau de retour [m]")
    ax4.set_title("Temps de retour vs Niveaux de retour de la nappe")
    ax4.legend(loc='lower right')
    ax4.grid(True, zorder=-4)

    # 5 . Histogram of max head difference per year
    ax5 = fig.add_subplot(326)
    years = list(max_head_diff_per_year.keys())
    max_head_diffs = list(max_head_diff_per_year.values())
    ax5.bar(years, max_head_diffs, color='royalblue', edgecolor='black', width=0.8, alpha=0.95, zorder=4)
    ax5.set_xlabel("Année")
    ax5.set_ylabel("Différence de niveau maximale [m]")
    ax5.set_title("Différence de niveau maximale par année et nombre de jours > Q95")
    ax5.grid(True, zorder=-4)

    # Label the number of days above Q95 for each year
    for year, max_diff in max_head_diff_per_year.items():
        days = days_above_Q95[year]
        ax5.text(year, max_diff, f'{days} j', ha='center', va='bottom', fontsize=10, color='black',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    # Add bars for specific years with the label "crue du Rhône"
    special_years = [2000, 2012, 2013, 2024]
    for year in special_years:
        ax5.bar(year, max(max_head_diffs), color='darkgreen', edgecolor='k', width=0.8, alpha=0.5, label="Crue du Rhône" if year == special_years[-1] else "")
        if year not in years:
            years.append(year)
    
    ax5.legend()
    ax5.set_xticks(years)
    ax5.set_xticklabels(years, rotation=45) 
    plt.tight_layout()

    if path_to_save is not False:
        plt.savefig(path_to_save, bbox_inches='tight')
    if to_plot:
        plt.show()
    else:
        plt.close()


def plot_max_serie_analysis_full(data_plot, name, x, y, TN, id, com, path_img, path_to_save=False, to_plot=True):
    data, extrem_event, u_gumbel, donnees_triees, slope, intercept, temps_retour, niveaux_retour, days_above_Q95, max_head_diff_per_year = data_plot
    
    fig = fig = plt.figure(figsize=(15, 18))   # A4 size in inches (21cm x 29.7cm)
    gs = gridspec.GridSpec(4, 2, figure=fig, height_ratios=[1, 1.2, 1, 1])  # Custom height ratios for better spacing

    # First row (2 subplots)
    info_ax = fig.add_subplot(gs[0, 0])
    info_ax.axis('off')  # Turn off the axis
    info_ax.text(0, 1, f"Nom : {name}", ha='left', va='top', fontsize=16, fontweight='bold')
    info_ax.text(0, 0.8, f"X : {x}", ha='left', va='top', fontsize=12)
    info_ax.text(0, 0.6, f"Y : {y}", ha='left', va='top', fontsize=12)
    info_ax.text(0, 0.4, f"TN : {TN}", ha='left', va='top', fontsize=12)    
    info_ax.text(0, 0.2, f"ID cantonal : {id}", ha='left', va='top', fontsize=12)
    info_ax.text(0, 0, f"Commune : {com}", ha='left', va='top', fontsize=12)

    # Map Axes (2/3 of the row)
    map_ax = fig.add_subplot(gs[0, 1:])
    map_ax.axis('off')  # Hide axis
    img = Image.open(path_img + 'map_{}.png'.format(name)).convert("RGB")
    map_ax.imshow(img, aspect='equal', extent=[0, img.size[0], 0, img.size[1]])
    #map_ax.imshow(img)
    print(path_img + 'map_{}.png'.format(name))

    # Add black border around the map
    rect = patches.Rectangle((0, 0), 1, 1, transform=map_ax.transAxes, 
                             linewidth=3, edgecolor='black', facecolor='none', zorder=10)
    map_ax.add_patch(rect)

    # Second row (1 full-width subplot)
    ax1 = fig.add_subplot(gs[1, :])
    ax1.plot(data.index, data['nvx_piezo'], label="Série Temporelle", color="k", alpha=0.6)
    ax1.scatter(extrem_event.index, extrem_event['nvx_piezo'], color="darkred", label="Max Annuel", zorder=3)
    ax1.set_title("Série temporelle piézométrique et maximums annuels piézo {}".format(name))
    ax1.set_xlabel("Date")
    ax1.set_ylabel("Niveau piézométrique")

    return_periods = np.array([10, 25, 50, 100])
    niveaux        = value_from_return_period(intercept, slope, return_periods)
    colors         = list(mcolors.LinearSegmentedColormap.from_list("", ["blue", "purple"])(np.linspace(0, 1, len(return_periods))))
    
    for rp, nvx, color in zip(return_periods, niveaux, colors):
        ax1.axhline(y=nvx, color=color, linestyle='--', linewidth=0.8, label=f'Temps de retour {rp} ans: {nvx:.2f} m')

    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    ax1.grid(True, zorder=-4)

    # Third row (2 subplots)
    ax2 = fig.add_subplot(gs[2, 0])
    ax2.scatter(u_gumbel, donnees_triees, label="Données", color='darkred')
    ax2.plot(u_gumbel, intercept + slope * u_gumbel, color="k", label=f"Régression: y = {slope:.2f}*u + {intercept:.2f}", linestyle='--')
    ax2.set_xlabel("Variable réduite u (Gumbel)")
    ax2.set_ylabel("Valeur Observée (m)")
    ax2.set_title("Ajustement des niveaux avec la loi de Gumbel")
    ax2.legend(loc='lower right')
    ax2.grid(True, zorder=-4)

    # 1. Distribution des valeurs par mois
    ax3 = fig.add_subplot(gs[2, 1])
    ax3.hist(extrem_event['Month'], bins=np.arange(1, 14) - 0.5, edgecolor='black', color='royalblue', alpha=0.95, width=0.8, zorder=4)
    ax3.set_xlabel("Mois")
    ax3.set_ylabel("Nombre d'événements")
    ax3.set_title("Distribution des maximums annuels par mois")
    ax3.set_xticks(np.arange(1, 13))
    ax3.set_xticklabels(['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Juin', 'Juil', 'Août', 'Sep', 'Oct', 'Nov', 'Déc'])
    ax3.grid(True, zorder=-4)

    # Fourth row (2 subplots)
    ax4 = fig.add_subplot(gs[3, 0])
    ax4.plot(temps_retour, niveaux_retour, label="Courbe de retour", color='k')
    important_events = extrem_event.sort_values('nvx_piezo')[-3:]
    ax4.scatter(return_period_from_value(intercept, slope, important_events['nvx_piezo'].values), important_events['nvx_piezo'], color='darkred', zorder=5, label="3 évènements les plus importants")

    # Save the 3 most important events and their return periods
    important_events['temps_retour'] = return_period_from_value(intercept, slope, important_events['nvx_piezo'].values)
    saved_important_events = important_events[['nvx_piezo', 'temps_retour']]

    # Annotate the scatter plot with the year of the events
    for i, event in important_events.iterrows():
        ax4.annotate(i.year, (return_period_from_value(intercept, slope, event['nvx_piezo']), event['nvx_piezo']), 
                     textcoords="offset points", xytext=(0,10), ha='center', 
                     bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    # 4 . Add horizontal line for TN (Terrain Naturel)
    ax4.hlines(y=TN, xmin=0, xmax=200, color='red', linestyle='--', label=f"Terrain naturel")
    ax4.set_xlabel("Temps de retour [années]")
    ax4.set_ylabel("Niveau de retour [m]")
    ax4.set_title("Temps de retour vs Niveaux de retour de la nappe")
    ax4.legend(loc='lower right')
    ax4.grid(True, zorder=-4)

    # 5 . Histogram of max head difference per year
    ax5 = fig.add_subplot(gs[3, 1])
    years          = list(max_head_diff_per_year.keys())
    max_head_diffs = list(max_head_diff_per_year.values())
    ax5.bar(years, max_head_diffs, color='royalblue', edgecolor='black', width=0.8, alpha=0.95, zorder=4, label="Max différence et nombre de jour > Q95")
    ax5.set_xlabel("Année")
    ax5.set_ylabel("Différence de niveau maximale [m]")
    ax5.set_title("Différence de niveau maximale par année et nombre de jours > Q95")
    ax5.grid(True, zorder=-4)

    # Label the number of days above Q95 for each year
    for year, max_diff in max_head_diff_per_year.items():
        days = days_above_Q95[year]
        ax5.text(year, max_diff, f'{days} j', ha='center', va='bottom', fontsize=10, color='black',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        
    # Add grey rectangle if the first data year is earlier than one of the special years
    special_years = np.array([2000, 2012, 2013, 2024])
    first_year    = min(years)
    
    if first_year > min(special_years):
        ax5.axvspan(min(special_years), first_year-1, facecolor='lightgrey', edgecolor='black', alpha=0.5, zorder=2, label='Pas de données piézométriques')
        
    for year in special_years:
        ax5.bar(year, max(max_head_diffs)+0.2*max(max_head_diffs), color='darkgreen', edgecolor='k', width=0.8, alpha=0.5, label="Crue du Rhône" if year == special_years[-1] else "")
        if year not in years:
            years.append(year)

    ax5.legend(loc='upper left')
    ax5.set_xticks(years)
    ax5.set_xticklabels(years, rotation=45)
    plt.tight_layout()
    
    if path_to_save:
        plt.savefig(path_to_save, bbox_inches='tight')
    if to_plot:
        plt.show()
        return saved_important_events
    else:
        plt.close()
        return saved_important_events
    
