import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks

def master_analysis():
    # Update paths if your script is not in the same folder as the .dat files
    files = {'Smooth': 'drag_smooth.dat', 'Bumpy': 'drag_bumpy.dat'}
    t_start = 15.0  # Skip startup transient (increased for Re=11,000)
    
    summary_data = []
    plt.figure(figsize=(12, 10))
    
    for i, (label, fname) in enumerate(files.items()):
        try:
            # Basilisk output format: t, Cd, Fp/ref, Fmu/ref
            # These are already normalized coefficients (Cdp and Cdv)
            data = pd.read_csv(fname, sep=' ', header=None, 
                             names=['t', 'Cd', 'Cdp', 'Cdv'], 
                             skipinitialspace=True)
        except FileNotFoundError:
            print(f"Error: {fname} not found. Skipping {label}...")
            continue

        # 1. Filter for Steady State
        steady = data[data['t'] > t_start].copy()
        
        # 2. Calculate Averages
        avg_cd = steady['Cd'].mean()
        avg_cdp = steady['Cdp'].mean()
        avg_cdv = steady['Cdv'].mean()
        
        # 3. Calculate Strouhal Number
        # Note: Drag oscillates at TWICE the shedding frequency (2f)
        peaks, _ = find_peaks(steady['Cd'], distance=10)
        if len(peaks) > 1:
            avg_period_drag = np.diff(steady['t'].iloc[peaks]).mean()
            freq_drag = 1.0 / avg_period_drag
            freq_shedding = freq_drag / 2.0  # The 2f correction
            st_num = (freq_shedding * 0.125) / 1.0  # St = (f*D)/U
        else:
            st_num = 0.0

        summary_data.append({
            'Case': label,
            'Total Cd': avg_cd,
            'Pressure Cd': avg_cdp,
            'Viscous Cd': avg_cdv,
            'Strouhal No.': st_num
        })

        # 4. Plotting Total Cd
        plt.subplot(2, 1, 1)
        plt.plot(steady['t'], steady['Cd'], label=f'{label} (Mean: {avg_cd:.3f})')
        plt.title('Total Drag Coefficient ($C_d$) at $Re=11,000$')
        plt.ylabel('$C_d$')
        plt.grid(True, alpha=0.3)
        plt.legend()

        # 5. Plotting Pressure vs Viscous (Stacked comparison)
        plt.subplot(2, 1, 2)
        plt.plot(steady['t'], steady['Cdp'], label=f'{label} Pressure component', alpha=0.7)
        plt.plot(steady['t'], steady['Cdv'], linestyle='--', label=f'{label} Viscous component')
        plt.title('Drag Component Breakdown')
        plt.xlabel('Time (non-dimensional)')
        plt.ylabel('Coefficient Value')
        plt.grid(True, alpha=0.3)
        plt.legend()

    # --- FINAL OUTPUTS ---
    plt.tight_layout()
    plt.savefig('btp_final_analysis.png', dpi=300)
    
    if summary_data:
        df_final = pd.DataFrame(summary_data)
        print("\n" + "="*70)
        print("                 FINAL BTP SUMMARY TABLE")
        print("="*70)
        print(df_final.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
        print("="*70)
        
        if len(summary_data) > 1:
            diff = ((summary_data[1]['Total Cd'] - summary_data[0]['Total Cd']) / summary_data[0]['Total Cd']) * 100
            print(f"\nResult: The Bumpy case has a {diff:+.2f}% change in total drag.")
            print("Note: Positive % means drag increased due to earlier separation.")
        print("="*70)
    
    plt.show()

if __name__ == "__main__":
    master_analysis()