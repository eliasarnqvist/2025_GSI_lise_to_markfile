# Elias Arnqvist, 2025-03-30

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import scipy.optimize as opt
import csv
import os
import re

#%% Settings

# Threshold for all yield peaks (in units of pps (parts-per-spill))
yield_threshold = 0.1
# Yield units conversion from pps to events per hour
yield_unit_conversion = 0.0014 * 3600

# Correct for too wide peaks?
wideness_correct_value = 14

# Markfile folder name
markfile_folder = "20250330_markfiles_N126_v2"
# Turn number in markfile
turn_number = -1
# Charge state of ions in markfile
ion_charge_state = -1

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

df_degraders = pd.read_csv("settings_degraders.csv", delimiter=';', dtype=str)


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
    
    for index, row in df_degraders.iterrows():
        this_setting = row["setting"]
        
        if this_setting == setting:
            this_degrader_number = row["degrader number"]
            this_degrader = row["degrader (mg/cm2)"]
            these_isotopes_main = row["main isotopes"].split(",")
            
            yields_main = []
            yields_overlapping = []
            
            for isotope, data in setting_data_sum.items():
                this_index = np.argmin(np.abs(data["x"] - float(this_degrader)))
                this_yield = data["y_corr"][this_index]
                
                if this_yield > yield_threshold:
                    if isotope in these_isotopes_main:
                        yields_main.append([isotope, this_yield])
                    else:
                        yields_overlapping.append([isotope, this_yield])
            
            setting_results[this_degrader_number] = {"yields_main" : yields_main,
                                                     "yields_overlapping" : yields_overlapping,
                                                     "degrader" : this_degrader}
            
    results[setting] = setting_results

#%% Save yields to markfiles

tab10 = plt.cm.tab10
hex_colors = ["0x" + clr.to_hex(tab10(i))[1:] for i in range(10)]

# Make folder structure
os.makedirs(markfile_folder, exist_ok=True)
for setting in unique_settings:
    os.makedirs(markfile_folder + "//" + str(setting), exist_ok=True)
    
    main_isotopes = []
    overlapping_isotopes = []
    
    for degrader_number, data in results[setting].items():
        degrader_value = data["degrader"]
        yields_main = data["yields_main"]
        yields_overlapping = data["yields_overlapping"]
        
        markfile_name_small = markfile_folder + "//" + str(setting) + "//degrader_" + degrader_number + "_" + degrader_value + ".mark.csv"
        with open(markfile_name_small, "w", newline="") as file_small:
            writer_small = csv.writer(file_small, delimiter=";")
            writer_small.writerow(["Name", "Turn", "Text", "Color"])
            
            for this_yield in yields_main:
                main_isotopes.append([this_yield[0], degrader_number])
                row = []
                isotope_name = re.sub(r"(\d{2,3})([A-Za-z]{1,2})", r"\2\1", this_yield[0])
                row.append("1" + isotope_name + ":-1e")
                row.append(str(turn_number))
                converted_yield = this_yield[1] * yield_unit_conversion
                row.append(this_yield[0] + "_" + f"{converted_yield:.2e}" + "_" + degrader_number)
                row.append(hex_colors[0])
                writer_small.writerow(row)
            
            for this_yield in yields_overlapping:
                overlapping_isotopes.append([this_yield[0], degrader_number])
                row = []
                isotope_name = re.sub(r"(\d{2,3})([A-Za-z]{1,2})", r"\2\1", this_yield[0])
                row.append("1" + isotope_name + ":-1e")
                row.append(str(turn_number))
                converted_yield = this_yield[1] * yield_unit_conversion
                row.append(this_yield[0] + "_" + f"{converted_yield:.2e}" + "_" + degrader_number)
                row.append(hex_colors[1])
                writer_small.writerow(row)

        markfile_name_all = markfile_folder + "//" + str(setting) + "//all_isotopes.mark.csv"
        with open(markfile_name_all, "w", newline="") as file_all:
            writer_all = csv.writer(file_all, delimiter=";")
            writer_all.writerow(["Name", "Turn", "Text", "Color"])
            
            for isotope in main_isotopes:
                row = []
                isotope_name = re.sub(r"(\d{2,3})([A-Za-z]{1,2})", r"\2\1", isotope[0])
                row.append("1" + isotope_name + ":-1e")
                row.append(str(turn_number))
                row.append(isotope[0] + "_" + isotope[1])
                row.append(hex_colors[0])
                writer_all.writerow(row)
            
            for isotope in overlapping_isotopes:
                row = []
                isotope_name = re.sub(r"(\d{2,3})([A-Za-z]{1,2})", r"\2\1", isotope[0])
                row.append("1" + isotope_name + ":-1e")
                row.append(str(turn_number))
                row.append(isotope[0] + "_" + isotope[1])
                row.append(hex_colors[1])
                writer_all.writerow(row)

#%% Save to a more human-readable file as well

info_file_name = markfile_folder + "//info.csv"
with open(info_file_name, "w", newline="") as file:
    writer = csv.writer(file, delimiter=";")
    writer.writerow(["setting", "degrader number", "degrader (mg/cm2)", "main isotopes (yield (/h))", "overlapping isotopes (yield (/h))"])
    
    for setting, data1 in results.items():
        for degrader_number, data2 in data1.items():
            row = []
            row.append(setting)
            row.append(degrader_number)
            row.append(data2["degrader"])
            isotopes_main = ""
            for isotope in data2["yields_main"]:
                isotopes_main += (isotope[0] + "(" + str(round(isotope[1], 2)) + ")_")
            isotopes_main = isotopes_main[:-1]
            row.append(isotopes_main)
            isotopes_overlapping = ""
            for isotope in data2["yields_overlapping"]:
                isotopes_overlapping += (isotope[0] + "(" + str(round(isotope[1], 2)) + ")_")
            isotopes_overlapping = isotopes_overlapping[:-1]
            row.append(isotopes_overlapping)
            
            writer.writerow(row)
