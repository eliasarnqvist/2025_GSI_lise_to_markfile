# 2025_GSI_lise_to_markfile

Written for the FRS Ion Catcher to produce markfiles from LISE++ data 

## Short instructions: 

1. In LISE++, go to Options > Plot options > Number of one-dimensional distributions. Change this number from 100 to at least 10000. 
2. Run LISE++ through Utilities > Range optimizer (Gas cell utility) > Plot already calculated isotopes... and wait for it to finish. When it has finished, click Save X-axis and plots to TXT-file. Place the file in the folder ```lise_files```. The name of the file must be the name of the setting you will use with the FRS. 
3. Manually fill in the FRS settings names and isotopes you want to measure in the file ```settings_isotopes.csv```. The setting name must match the corresponding LISE++ file from above. 
4. Run ```get_data_and_plot.py```. This should produce plots of the isotopes to measure, for each setting. Both LISE++ data and corrected data is shown. The correction was done by fitting Gaussians to peaks and forcing them to not bee too wide. 
5. Analyze the plots from above and decide what degrader thickness you want to use. Manually add this information to ```settings_degraders.csv```. Here you must also list all isotopes you want to measure for a given degrader setting. 
6. Open ```get_data_and_make_markfile.py``` and change settings (folder name, turn number, ion charge, etc.) as you wish. Then run ```get_data_and_make_markfile.py```. This produces a folder with a name of your choosing containing all markfiles. 
7. Load these markfiles as you wish in TOFControl. 

## Limitations: 

- LISE++ sometimes makes too wide distributions. This is apparently a bug. Be careful. The way this script corrects for this is not perfect. 
- You still need to manually inspect what degrader setting is best by looking at the plot. This is difficult to make automatic, as the LISE++ data sometimes does not match the correction. 
- The code is not written with efficiency in mind, but this is not really an issue. 

For the FRS Ion Catcher

Elias Arnqvist 2025