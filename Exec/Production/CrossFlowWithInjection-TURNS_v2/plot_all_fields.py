#!/usr/bin/env python3
import yt
import os
import glob
import sys

def create_output_folder(script_dir):
    """Create output folder in the script's directory"""
    output_folder = os.path.join(script_dir, "output")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    return output_folder

def plot_all_fields(plt_path, output_folder):
    """Plot all available fields from a plotfile"""
    ds = yt.load(plt_path)
    plt_num = os.path.basename(plt_path)
    
    field_list = [f for f in ds.field_list if f[0] == 'boxlib']
    
    for field in field_list:
        field_name = field[1]
        output_path = os.path.join(output_folder, f"{plt_num}_{field_name}.png")
        
        try:
            slc = yt.SlicePlot(ds, "z", field)
            
            # Apply field-specific settings
            if "temp" in field_name.lower():
                slc.set_cmap(field, "hot")
                slc.set_log(field, False)
            elif "vel" in field_name.lower():
                slc.set_cmap(field, "viridis")
            elif "pres" in field_name.lower():
                slc.set_cmap(field, "plasma")
                slc.set_log(field, True)
            else:
                slc.set_cmap(field, "jet")
            
            # Add EB visualization if available
            if ("boxlib", "vfrac") in ds.field_list:
                slc.annotate_contour(("boxlib", "vfrac"), levels=[0.5], clim=(0,1),
                                   plot_args={"colors": "white", "linewidths": 1})
            
            slc.annotate_title(f"{field_name}")
            slc.save(output_path)
            print(f"Saved {field_name} plot to {output_path}")
            
        except Exception as e:
            print(f"Failed to plot {field_name}: {str(e)}")

def process_all_plt_files(script_dir):
    """Process all plt files in the script's directory"""
    output_folder = create_output_folder(script_dir)
    plt_files = sorted(glob.glob(os.path.join(script_dir, 'plt[0-9][0-9][0-9][0-9][0-9]')))
    
    if not plt_files:
        print(f"No plt files found in {script_dir}")
        return
    
    for plt_file in plt_files:
        print(f"\nProcessing {os.path.basename(plt_file)}...")
        plot_all_fields(plt_file, output_folder)

if __name__ == "__main__":
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    print(f"Looking for PLT files in: {script_dir}")
    process_all_plt_files(script_dir)
    print("\nAll plots saved to 'output' folder in the script directory.")