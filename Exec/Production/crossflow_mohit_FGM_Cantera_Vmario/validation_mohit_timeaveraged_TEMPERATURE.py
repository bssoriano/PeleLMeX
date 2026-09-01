import yt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import re

# --- PARAMETERS ---
D = 0.0762  # Diameter in meters
first_pf = 'plt01011'
last_pf = 'plt30000'

####################### PLOT 1D ##############################################################
csv_file = 'mohit_fig4aT.csv'
# 1. Automatically generate the list of plotfiles in the range
def get_plotfiles_in_range(start_name, end_name):
    # Extract numbers from the names (e.g., 'plt00203' -> 203)
    start_num = int(re.search(r'\d+', start_name).group())
    end_num = int(re.search(r'\d+', end_name).group())
    
    # List all directories starting with 'plt'
    all_files = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('plt')]
    
    selected = []
    for f in all_files:
        match = re.search(r'\d+', f)
        if match:
            num = int(match.group())
            if start_num <= num <= end_num:
                selected.append(f)
    
    return sorted(selected)

plotfiles = get_plotfiles_in_range(first_pf, last_pf)
print(f"Found {len(plotfiles)} plotfiles in range: {plotfiles}")
#plotfiles = ['plt01131', 'plt02011', 'plt03001','plt04374','plt06040']

# Lineout parameters
target_y = 10.5*D + 1.0*D + 0.001 # Fixed axial height in meters
target_z = 0.0     # Fixed z coordinate

# 1. Load the Reference CSV
# Assuming Column 0 = Temperature, Column 1 = x/D
df_ref = pd.read_csv(csv_file, sep=';', decimal=',', header=None)
df_ref.columns = ['Temperature', 'x_D']

# Define a common radial grid for interpolation (based on CSV x/D range)
x_common_norm = np.linspace(df_ref['x_D'].min(), df_ref['x_D'].max(), 500)
x_common_phys = x_common_norm * D

# 2. Extract and Process PeleLMex Data along the X-axis
all_temps_interp = []

print(f"Extracting radial profiles at y = {target_y} m...")
for pf in plotfiles:
    try:
        ds = yt.load(pf)
        
        # ortho_ray(axis, coords) 
        # axis 0 = X-axis. coords = (fixed_y, fixed_z)
        ray = ds.ortho_ray(0, (target_y, target_z))
        
        # Sort and clean data
        raw_x = np.array(ray["x"])
        raw_t = np.array(ray["temp"])
        sort_idx = np.argsort(raw_x)
        raw_x, raw_t = raw_x[sort_idx], raw_t[sort_idx]
        
        # Remove duplicates
        raw_x, unique_idx = np.unique(raw_x, return_index=True)
        raw_t = raw_t[unique_idx]
        
        # Interpolate onto the common grid
        f_interp = interp1d(raw_x, raw_t, bounds_error=False, fill_value="extrapolate")
        all_temps_interp.append(f_interp(x_common_phys))
        
        print(f" - {pf} processed.")
    except Exception as e:
        print(f" - Error in {pf}: {e}")

# 3. Calculate Mean and RMS
mean_temp = np.mean(all_temps_interp, axis=0)
std_temp = np.std(all_temps_interp, axis=0)

# 4. Plotting
plt.figure(figsize=(8, 10))

# Experimental/Reference Data (Points)
plt.plot(df_ref['Temperature'], df_ref['x_D'], 
            color='blue',linestyle= "--", alpha=0.6, label='Reference Fig. 4 (Mohit, 2025)')

# Mean Simulation Profile (Line)
plt.plot(mean_temp, x_common_norm, 
         color='green', linewidth=2, label='PeleLMex Time-averaged')

# RMS Fluctuations (Shaded Area)
plt.fill_betweenx(x_common_norm, mean_temp - std_temp, mean_temp + std_temp, 
                  color='green', alpha=0.2, label='PeleLMex RMS')

# Formatting
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('x / D', fontsize=12)
#plt.title(f'Temperature Profile at y = {target_y} m\n(Time-Average across {len(plotfiles)} snapshots)')
#plt.grid(True, linestyle=':', alpha=0.7)
plt.xlim(0,2000)
plt.ylim(0,15)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('mohit_fig4aT.png', dpi=300)
#plt.show()
###################################################################################################
#
####################### PLOT 6D ##############################################################
csv_file = 'mohit_fig4bT.csv'
# 1. Automatically generate the list of plotfiles in the range
def get_plotfiles_in_range(start_name, end_name):
    # Extract numbers from the names (e.g., 'plt00203' -> 203)
    start_num = int(re.search(r'\d+', start_name).group())
    end_num = int(re.search(r'\d+', end_name).group())
    
    # List all directories starting with 'plt'
    all_files = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('plt')]
    
    selected = []
    for f in all_files:
        match = re.search(r'\d+', f)
        if match:
            num = int(match.group())
            if start_num <= num <= end_num:
                selected.append(f)
    
    return sorted(selected)

plotfiles = get_plotfiles_in_range(first_pf, last_pf)
print(f"Found {len(plotfiles)} plotfiles in range: {plotfiles}")
#plotfiles = ['plt01131', 'plt02011', 'plt03001','plt04374','plt06040']

# Lineout parameters
target_y = 10.0*D + 6.0*D  # Fixed axial height in meters
target_z = 0.0     # Fixed z coordinate

# 1. Load the Reference CSV
# Assuming Column 0 = Temperature, Column 1 = x/D
df_ref = pd.read_csv(csv_file, sep=';', decimal=',', header=None)
df_ref.columns = ['Temperature', 'x_D']

# Define a common radial grid for interpolation (based on CSV x/D range)
x_common_norm = np.linspace(df_ref['x_D'].min(), df_ref['x_D'].max(), 500)
x_common_phys = x_common_norm * D

# 2. Extract and Process PeleLMex Data along the X-axis
all_temps_interp = []

print(f"Extracting radial profiles at y = {target_y} m...")
for pf in plotfiles:
    try:
        ds = yt.load(pf)
        
        # ortho_ray(axis, coords) 
        # axis 0 = X-axis. coords = (fixed_y, fixed_z)
        ray = ds.ortho_ray(0, (target_y, target_z))
        
        # Sort and clean data
        raw_x = np.array(ray["x"])
        raw_t = np.array(ray["temp"])
        sort_idx = np.argsort(raw_x)
        raw_x, raw_t = raw_x[sort_idx], raw_t[sort_idx]
        
        # Remove duplicates
        raw_x, unique_idx = np.unique(raw_x, return_index=True)
        raw_t = raw_t[unique_idx]
        
        # Interpolate onto the common grid
        f_interp = interp1d(raw_x, raw_t, bounds_error=False, fill_value="extrapolate")
        all_temps_interp.append(f_interp(x_common_phys))
        
        print(f" - {pf} processed.")
    except Exception as e:
        print(f" - Error in {pf}: {e}")

# 3. Calculate Mean and RMS
mean_temp = np.mean(all_temps_interp, axis=0)
std_temp = np.std(all_temps_interp, axis=0)

# 4. Plotting
plt.figure(figsize=(8, 10))

# Experimental/Reference Data (Points)
plt.plot(df_ref['Temperature'], df_ref['x_D'], 
            color='blue',linestyle= "--", alpha=0.6, label='Reference Fig. 4 (Mohit, 2025)')

# Mean Simulation Profile (Line)
plt.plot(mean_temp, x_common_norm, 
         color='green', linewidth=2, label='PeleLMex Time-averaged')

# RMS Fluctuations (Shaded Area)
plt.fill_betweenx(x_common_norm, mean_temp - std_temp, mean_temp + std_temp, 
                  color='green', alpha=0.2, label='PeleLMex RMS')

# Formatting
plt.xlabel('Temperature [K]', fontsize=12)
plt.ylabel('x / D', fontsize=12)
#plt.title(f'Temperature Profile at y = {target_y} m\n(Time-Average across {len(plotfiles)} snapshots)')
#plt.grid(True, linestyle=':', alpha=0.7)
plt.xlim(0,2000)
plt.ylim(0,15)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('mohit_fig4bT.png', dpi=300)
plt.show()
###################################################################################################



####################### PLOT 11D ##############################################################
#csv_file = 'mohit_fig4cT.csv'
## 1. Automatically generate the list of plotfiles in the range
#def get_plotfiles_in_range(start_name, end_name):
#    # Extract numbers from the names (e.g., 'plt00203' -> 203)
#    start_num = int(re.search(r'\d+', start_name).group())
#    end_num = int(re.search(r'\d+', end_name).group())
#    
#    # List all directories starting with 'plt'
#    all_files = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('plt')]
#    
#    selected = []
#    for f in all_files:
#        match = re.search(r'\d+', f)
#        if match:
#            num = int(match.group())
#            if start_num <= num <= end_num:
#                selected.append(f)
#    
#    return sorted(selected)
#
#plotfiles = get_plotfiles_in_range(first_pf, last_pf)
#print(f"Found {len(plotfiles)} plotfiles in range: {plotfiles}")
##plotfiles = ['plt01131', 'plt02011', 'plt03001','plt04374','plt06040']
#
## Lineout parameters
#target_y = 10.0*D + 11.0*D  # Fixed axial height in meters
#target_z = 0.0     # Fixed z coordinate
#
## 1. Load the Reference CSV
## Assuming Column 0 = Temperature, Column 1 = x/D
#df_ref = pd.read_csv(csv_file, sep=';', decimal=',', header=None)
#df_ref.columns = ['Temperature', 'x_D']
#
## Define a common radial grid for interpolation (based on CSV x/D range)
#x_common_norm = np.linspace(df_ref['x_D'].min(), df_ref['x_D'].max(), 500)
#x_common_phys = x_common_norm * D
#
## 2. Extract and Process PeleLMex Data along the X-axis
#all_temps_interp = []
#
#print(f"Extracting radial profiles at y = {target_y} m...")
#for pf in plotfiles:
#    try:
#        ds = yt.load(pf)
#        
#        # ortho_ray(axis, coords) 
#        # axis 0 = X-axis. coords = (fixed_y, fixed_z)
#        ray = ds.ortho_ray(0, (target_y, target_z))
#        
#        # Sort and clean data
#        raw_x = np.array(ray["x"])
#        raw_t = np.array(ray["temp"])
#        sort_idx = np.argsort(raw_x)
#        raw_x, raw_t = raw_x[sort_idx], raw_t[sort_idx]
#        
#        # Remove duplicates
#        raw_x, unique_idx = np.unique(raw_x, return_index=True)
#        raw_t = raw_t[unique_idx]
#        
#        # Interpolate onto the common grid
#        f_interp = interp1d(raw_x, raw_t, bounds_error=False, fill_value="extrapolate")
#        all_temps_interp.append(f_interp(x_common_phys))
#        
#        print(f" - {pf} processed.")
#    except Exception as e:
#        print(f" - Error in {pf}: {e}")
#
## 3. Calculate Mean and RMS
#mean_temp = np.mean(all_temps_interp, axis=0)
#std_temp = np.std(all_temps_interp, axis=0)
#
## 4. Plotting
#plt.figure(figsize=(8, 10))
#
## Experimental/Reference Data (Points)
#plt.plot(df_ref['Temperature'], df_ref['x_D'], 
#            color='blue',linestyle= "--", alpha=0.6, label='Reference Fig. 4 (Mohit, 2025)')
#
## Mean Simulation Profile (Line)
#plt.plot(mean_temp, x_common_norm, 
#         color='green', linewidth=2, label='PeleLMex Time-averaged')
#
## RMS Fluctuations (Shaded Area)
#plt.fill_betweenx(x_common_norm, mean_temp - std_temp, mean_temp + std_temp, 
#                  color='green', alpha=0.2, label='PeleLMex RMS')
#
## Formatting
#plt.xlabel('Temperature [K]', fontsize=12)
#plt.ylabel('x / D', fontsize=12)
##plt.title(f'Temperature Profile at y = {target_y} m\n(Time-Average across {len(plotfiles)} snapshots)')
##plt.grid(True, linestyle=':', alpha=0.7)
#plt.xlim(0,2000)
#plt.ylim(0,15)
#plt.legend(loc='best')
#
#plt.tight_layout()
#plt.savefig('mohit_fig4cT.png', dpi=300)
#plt.show()
###################################################################################################
