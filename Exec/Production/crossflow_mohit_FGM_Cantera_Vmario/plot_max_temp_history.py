#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np
import os
import sys 
import glob

# Suppress verbose yt loading logs to keep the terminal output clean
yt.funcs.mylog.setLevel(40)

def extract_max_temp_history(base_dir):
    """
    Loops through all plt* folders, extracts simulation time and max temperature
    strictly on the Z-center slice.
    """
    plot_pattern = os.path.join(base_dir, "plt*")
    plot_dirs = sorted([f for f in glob.glob(plot_pattern) if os.path.isdir(f)])

    if not plot_dirs:
        print(f"Error: No 'plt*' directories found in {base_dir}")
        sys.exit(1)

    print(f"--- Found {len(plot_dirs)} plotfiles to process ---")

    times = []
    max_temps = []
    field_name = ('boxlib', 'temp')

    for plt_path in plot_dirs:
        plot_name = os.path.basename(plt_path)
        try:
            ds = yt.load(plt_path)
            
            if field_name not in ds.field_list:
                print(f"  [WARNING] Field '{field_name}' not found in {plot_name}. Skipping.")
                continue
            
            # Extract current simulation time (seconds)
            sim_time = float(ds.current_time.v)
            
            # --- MODIFIED SECTION ---
            # Get the exact numerical coordinate for the center of the Z-axis
            z_center = ds.domain_center[2] 
            
            # Create a 2D data object cutting through that exact Z coordinate
            slc = ds.slice('z', z_center)
            
            # Find the maximum value strictly on this 2D slice
            max_t = float(slc[field_name].max())
            # ------------------------
            
            times.append(sim_time)
            max_temps.append(max_t)
            
            print(f"Processed {plot_name}: Time = {sim_time:.6e} s | Z-Slice Max Temp = {max_t:.2f} K")
            
        except Exception as e:
            print(f"  [ERROR] Failed to process {plot_name}: {e}")

    return np.array(times), np.array(max_temps)

def plot_history(times, max_temps, save_path):
    """
    Generates and saves the Time vs Max Temperature plot.
    """
    if len(times) == 0:
        print("No valid data collected. Plot creation aborted.")
        return

    plt.figure(figsize=(9, 5.5))
    
    # Plotting line with small indicators for each data point
    plt.plot(times, max_temps, marker='o', linestyle='-', color='darkorange', linewidth=2, markersize=4, label='Z-Slice Max')
    
    # Formatting
    plt.title("Maximum Z-Slice Temperature Evolution", fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Simulation Time (s)", fontsize=12)
    plt.ylabel("Maximum Temperature (K)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    # Use scientific notation for time axis if steps are very tiny
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot to: {save_path}")
    plt.show()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract and plot Max Z-Slice Temperature over time from AMReX plotfiles.')
    parser.add_argument("base_dir", help="Path to the directory containing your pltXXXXX folders")
    parser.add_argument("-o", "--output", default="max_temp_vs_time_zslice.png", help="Output filename for the image")
    args = parser.parse_args()

    # 1. Extract data
    times, max_temps = extract_max_temp_history(args.base_dir)
    
    # 2. Save raw data to a text file
    data_file = os.path.join(args.base_dir, "max_temp_history_zslice.txt")
    np.savetxt(data_file, np.column_stack((times, max_temps)), header="Time(s) MaxTemp_ZSlice(K)", fmt="%.6e")
    print(f"\nSaved raw data table to: {data_file}")

    # 3. Generate the final plot
    output_plot_path = os.path.join(args.base_dir, args.output)
    plot_history(times, max_temps, output_plot_path)
    
    print("\n--- Execution Complete ---")