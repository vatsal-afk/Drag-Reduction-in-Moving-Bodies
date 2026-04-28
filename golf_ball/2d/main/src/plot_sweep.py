import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

# The directory where run_range.sh saves the simulation results
BASE_DIR = "../results/parametric_sweep"

def main():
    print(f"Scanning for data in {BASE_DIR}...")
    
    # Use glob to find all the created 'drag_bumpy_*.dat' files
    drag_files = glob.glob(f"{BASE_DIR}/Re_*/drag_*.dat")
    
    if not drag_files:
        print(f"No drag files found! Did the simulation run properly and save inside {BASE_DIR}?")
        return
        
    mean_cd_dict = {}

    # Setup the first plot for time history
    plt.figure(figsize=(10, 6))

    for fpath in drag_files:
        filename = os.path.basename(fpath)
        
        # Regex to extract the Reynolds number from e.g. "drag_bumpy_Re5000.dat"
        match = re.search(r'Re(\d+(\.\d+)?)', filename)
        if not match:
            continue
            
        Re = float(match.group(1))
        
        try:
            # Data columns are: Time | Cd | Fp (Pressure Drag) | Fmu (Viscous Drag)
            data = np.loadtxt(fpath)
            
            # If the file only has one row, reshape it to 2D
            if data.ndim == 1:
                data = data.reshape(1, -1)
                
            t_vals = data[:, 0]
            cd_vals = data[:, 1]
            
            # Plot Cd vs Time 
            plt.plot(t_vals, cd_vals, label=f"Re = {int(Re) if Re.is_integer() else Re}")
            
            # To get a stable "Mean Cd", we average only the second half of the 
            # simulation data to ignore the initial transient spike when flow starts
            start_steady = len(cd_vals) // 2
            mean_cd = np.mean(cd_vals[start_steady:]) 
            
            mean_cd_dict[Re] = mean_cd
            
        except Exception as e:
            print(f"Error processing {fpath}: {e}")

    # Plot 1: Cd vs Time for all Re
    plt.xlabel("Time")
    plt.ylabel("Drag Coefficient ($C_d$)")
    plt.title("Transient Drag Coefficient vs Time")
    plt.legend()
    plt.grid(True)
    
    plot1_path = os.path.join(BASE_DIR, "plot_Cd_vs_Time.png")
    plt.savefig(plot1_path, dpi=300)
    print(f"Saved: {plot1_path}")
    plt.close()

    # Plot 2: Mean Cd vs  Reynolds Number
    if mean_cd_dict:
        # Sort the dictionary keys (the Reynolds numbers) in ascending order
        sorted_res = sorted(mean_cd_dict.keys())
        sorted_cds = [mean_cd_dict[r] for r in sorted_res]
        
        plt.figure(figsize=(8, 5))
        plt.plot(sorted_res, sorted_cds, marker='o', linestyle='-', color='b')
        plt.xscale('log') # Use log scale for Re since it goes from 100 to 50000+
        plt.xlabel("Reynolds Number (Re)")
        plt.ylabel("Mean Drag Coefficient ($C_d$)")
        plt.title("Mean Drag Coefficient vs Reynolds Number")
        plt.grid(True, which="both", ls="--", alpha=0.5)
        
        plot2_path = os.path.join(BASE_DIR, "plot_Mean_Cd_vs_Re.png")
        plt.savefig(plot2_path, dpi=300)
        print(f"Saved: {plot2_path}")
        plt.close()

if __name__ == "__main__":
    main()
