# Elias Arnqvist, 2025-03-30

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

#%% Settings

# Correct for too wide peaks?
wideness_correct_value = 14

#%% Functions

def gaussian(x, a, b, c):
    return a * np.exp(-((x - b) ** 2) / (2 * c ** 2))

def curve_correction (x, y):
    if np.sum(y) > 0:
        mean = np.average(x, weights=y)
        std = np.sqrt(np.average((x - mean) ** 2, weights=y))
        initial_guess = [np.max(y), mean, std]
        popt, pcov = opt.curve_fit(gaussian, x, y, p0=initial_guess)
        
        # Now correct for the width
        popt[0] = popt[0] * popt[2] / wideness_correct_value
        popt[2] = wideness_correct_value
        
        return gaussian(x, *popt)
    else:
        return y

#%% Import names of isotopes and FRS settings

# Import the names of settings we use and isotopes we measure
df_isotopes = pd.read_csv("settings_isotopes.csv", delimiter=';', dtype=str)
unique_settings = df_isotopes['setting'].unique().tolist()

#%% Import the data from LISE++

data_lise = {}
for setting in unique_settings:
    file_name = "lise_files\\" + setting + ".txt"
    df_file = pd.read_csv(file_name, delimiter='\t', skiprows=4)
    
    data_lise[setting] = df_file

isotope_data = {}
for setting in unique_settings:
    setting_data = {}
    
    x = data_lise[setting].iloc[:, 0].to_numpy()
    for yield_curve in data_lise[setting].iloc[:, 1::2].columns:
        isotope_name = yield_curve.split()[0]
        y = data_lise[setting][yield_curve].to_numpy()
        y_corr = curve_correction(x, y)
        
        if isotope_name not in setting_data.keys():
            setting_data[isotope_name] = {"x" : [x],
                                          "y" : [y],
                                          "y_corr" : [y_corr]}
        else:
            setting_data[isotope_name]["x"].append(x)
            setting_data[isotope_name]["y"].append(y)
            setting_data[isotope_name]["y_corr"].append(y_corr)
    isotope_data[setting] = setting_data

isotope_data_sum = {}
for setting in unique_settings:
    setting_data = isotope_data[setting]
    setting_data_sum = {}
    
    for isotope, data in setting_data.items():
        x = data["x"][0]
        y = np.sum(data["y"], axis=0)
        y_corr = np.sum(data["y_corr"], axis=0)
        
        setting_data_sum[isotope] = {"x" : x,
                                     "y" : y,
                                     "y_corr" : y_corr}
    isotope_data_sum[setting] = setting_data_sum

results = {}
for setting in unique_settings:
    setting_data_sum = isotope_data_sum[setting]
    
    setting_results = {}
    
    for index, row in df_isotopes.iterrows():
        this_setting = row["setting"]
        this_isotope = row["isotope"]
        
        if this_setting == setting:
            for isotope, data in setting_data_sum.items():
                if isotope == this_isotope:
                    y_max = np.max(data["y"])
                    x_at_y_max = data["x"][np.argmax(data["y"])]
                    y_corr_max = np.max(data["y_corr"])
                    x_at_y_corr_max = data["x"][np.argmax(data["y_corr"])]
                    
                    setting_results[isotope] = {"y_max" : [x_at_y_max, y_max],
                                                "y_corr_max" : [x_at_y_corr_max, y_corr_max]}
                
    results[setting] = setting_results

#%% Plot things so we can check if it is accurate

plt.close('all')
inch_to_mm = 25.4
# color = plt.cm.tab10
color = plt.cm.tab20

for setting in unique_settings:
    fig, ax = plt.subplots(figsize=(160/inch_to_mm, 100/inch_to_mm))
    i = 0
    
    setting_results = results[setting]
    for isotope, data in setting_results.items():
        y_max = data ["y_max"]
        y_corr_max = data["y_corr_max"]
        
        setting_data_sum = isotope_data_sum[setting][isotope]
        x = setting_data_sum["x"]
        y = setting_data_sum["y"]
        y_corr = setting_data_sum["y_corr"]
        
        ax.plot(x, y, ls="--", color=color(i))
        ax.plot(x, y_corr, ls="-", color=color(i), label=isotope + " (" + str(y_max[0]) + ", corr: " + str(y_corr_max[0]) + ")")
        
        ax.plot([y_max[0], y_max[0]], [0, y_max[1]], ls="--", marker="v", color=color(i))
        ax.plot([y_corr_max[0], y_corr_max[0]], [0, y_corr_max[1]], ls="-", marker="^", color=color(i))
        
        i += 1
    
    ax.set_title(setting)
    ax.legend(frameon=False, ncols=2)
    plt.tight_layout(pad = 0.2)
