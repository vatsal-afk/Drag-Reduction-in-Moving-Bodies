import matplotlib.pyplot as plt
import pandas as pd
import os

def plot_zoomed_results():
    smooth_file = 'drag_smooth.dat'
    bumpy_file = 'drag_bumpy.dat'

    if not os.path.exists(smooth_file) or not os.path.exists(bumpy_file):
        print("Error: Files not found.")
        return

    cols = ['time', 'Cd', 'Fp', 'Fmu']
    smooth = pd.read_csv(smooth_file, sep=' ', header=None, names=cols)
    bumpy = pd.read_csv(bumpy_file, sep=' ', header=None, names=cols)

    # Filter out the initial startup spike
    t_start = 10.0
    s_data = smooth[smooth['time'] > t_start]
    b_data = bumpy[bumpy['time'] > t_start]

    plt.figure(figsize=(12, 6))

    # Plotting with distinct colors and slightly thicker lines
    plt.plot(s_data['time'], s_data['Cd'], label='Smooth Cylinder', color='#1f77b4', linewidth=2)
    plt.plot(b_data['time'], b_data['Cd'], label='Bumpy (Dimpled)', color='#ff7f0e', linewidth=2)

    # --- THE FIX: TIGHT SCALING ---
    # We find the min and max of the ACTUAL data and add a 5% margin
    all_values = pd.concat([s_data['Cd'], b_data['Cd']])
    data_min = all_values.min()
    data_max = all_values.max()
    padding = (data_max - data_min) * 0.5 if data_max != data_min else 0.05
    
    # Set the Y-axis to zoom in specifically on the results
    plt.ylim(data_min - padding, data_max + padding)

    plt.title('Zoomed View: Drag Coefficient Comparison', fontsize=14)
    plt.xlabel('Time (t)', fontsize=12)
    plt.ylabel('Drag Coefficient (Cd)', fontsize=12)
    plt.legend(loc='best')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('zoomed_drag_plot.png', dpi=300)
    print(f"Mean Smooth: {s_data['Cd'].mean():.4f}")
    print(f"Mean Bumpy:  {b_data['Cd'].mean():.4f}")
    plt.show()

if __name__ == "__main__":
    plot_zoomed_results()