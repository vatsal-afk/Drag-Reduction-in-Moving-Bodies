"""
Step 4: Drag Comparison — Golf Ball vs Smooth Sphere
=====================================================
Reads drag_smooth.dat and drag_golf.dat produced by the Basilisk simulations
and generates:
  - Time-series plot of Cd for both cases
  - Pressure vs viscous drag breakdown
  - Drag reduction statistics
  - drag_comparison.png

Usage:
    pip install numpy matplotlib scipy
    python step4_analyze.py

Expects in the same directory:
    drag_smooth.dat    (from smooth_sphere simulation)
    drag_golf.dat      (from golf_ball simulation)
    result_smooth.txt  (summary file from smooth_sphere)
    result_golf.txt    (summary file from golf_ball)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import sys


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_drag(filename):
    """Load drag data file into dict of arrays."""
    data = np.loadtxt(filename, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {
        't'            : data[:, 0],
        'Cd'           : data[:, 1],
        'Cl'           : data[:, 2],
        'Cd_pressure'  : data[:, 3],
        'Cd_viscous'   : data[:, 4],
    }


def moving_average(x, window=20):
    """Simple moving average for smoothing noisy Cd time series."""
    return np.convolve(x, np.ones(window) / window, mode='same')


def read_summary(filename):
    """Read Re and mean Cd from result_*.txt summary."""
    result = {}
    with open(filename) as f:
        for line in f:
            k, v = line.strip().split('=')
            result[k] = float(v)
    return result


def compute_stats(data, t_stat=40.):
    """Compute mean, std, min, max of Cd for t >= t_stat."""
    mask = data['t'] >= t_stat
    cd = data['Cd'][mask]
    if len(cd) == 0:
        return None
    return {
        'mean'  : np.mean(cd),
        'std'   : np.std(cd),
        'min'   : np.min(cd),
        'max'   : np.max(cd),
        'n'     : len(cd),
        'Cd_p'  : np.mean(data['Cd_pressure'][mask]),
        'Cd_v'  : np.mean(data['Cd_viscous'][mask]),
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    t_stat = 40.0   # start averaging after this time

    # ── Load data ──────────────────────────────────────────────────────────────
    files_needed = ['drag_smooth.dat', 'drag_golf.dat']
    missing = [f for f in files_needed if not Path(f).exists()]
    if missing:
        print(f"ERROR: Missing files: {missing}")
        print("Run smooth_sphere and golf_ball simulations first.")
        sys.exit(1)

    smooth = load_drag('drag_smooth.dat')
    golf   = load_drag('drag_golf.dat')

    # Read Reynolds number from summary (or fall back to inferring from data)
    Re = 300.
    if Path('result_smooth.txt').exists():
        s = read_summary('result_smooth.txt')
        Re = s['Re']

    # ── Compute statistics ─────────────────────────────────────────────────────
    stats_s = compute_stats(smooth, t_stat)
    stats_g = compute_stats(golf,   t_stat)

    if stats_s is None or stats_g is None:
        print(f"WARNING: Not enough data after t={t_stat}. "
              f"Reduce t_stat or run longer simulation.")
        t_stat = min(smooth['t'].max(), golf['t'].max()) * 0.3
        stats_s = compute_stats(smooth, t_stat)
        stats_g = compute_stats(golf,   t_stat)

    Cd_s = stats_s['mean']
    Cd_g = stats_g['mean']
    drag_reduction = 100. * (Cd_s - Cd_g) / Cd_s

    # ── Print report ───────────────────────────────────────────────────────────
    print("\n" + "="*55)
    print(f"  DRAG COMPARISON: Golf Ball vs Smooth Sphere")
    print(f"  Reynolds number: Re = {Re:.0f}")
    print("="*55)
    print(f"\n{'':25s} {'Smooth Sphere':>16s}  {'Golf Ball':>12s}")
    print(f"  {'─'*51}")
    print(f"  {'Cd (mean)':24s} {Cd_s:>16.4f}  {Cd_g:>12.4f}")
    print(f"  {'Cd (std dev)':24s} {stats_s['std']:>16.4f}  {stats_g['std']:>12.4f}")
    print(f"  {'Cd_pressure':24s} {stats_s['Cd_p']:>16.4f}  {stats_g['Cd_p']:>12.4f}")
    print(f"  {'Cd_viscous':24s} {stats_s['Cd_v']:>16.4f}  {stats_g['Cd_v']:>12.4f}")
    print(f"  {'Samples (t≥{:.0f})'.format(t_stat):24s} {stats_s['n']:>16d}  {stats_g['n']:>12d}")
    print(f"\n  Drag reduction (golf vs smooth) : {drag_reduction:+.1f}%")
    if drag_reduction > 0:
        print(f"  → Golf ball has LESS drag at Re={Re:.0f}")
        print(f"    (dimples trip boundary layer, delay separation)")
    elif drag_reduction < 0:
        print(f"  → Golf ball has MORE drag at Re={Re:.0f}")
        print(f"    (below critical Re, dimples increase drag — expected!)")
        print(f"    Critical Re for drag crossover ≈ 4×10⁴")
    print()

    # ── Physical interpretation ─────────────────────────────────────────────────
    print("  Physical context:")
    print(f"  At Re={Re:.0f}:")
    if Re < 4e4:
        print("    • Flow is laminar around a smooth sphere")
        print("    • Dimples force premature boundary layer transition")
        print("    • This INCREASES drag below the critical Re")
        print("    • Real golf balls travel at Re~1×10⁵ where dimples help")
    else:
        print("    • Flow is turbulent / transitional around smooth sphere")
        print("    • Dimples lower separation point → smaller wake → less drag")
        print("    • This is the 'drag crisis' regime golf balls exploit")
    print()

    # ── Literature comparison ───────────────────────────────────────────────────
    lit = {300: (0.65, None), 1000: (0.47, None),
           3700: (0.40, None), 100000: (0.47, 0.27)}
    if int(Re) in lit:
        cd_smooth_lit, cd_golf_lit = lit[int(Re)]
        print(f"  Literature values at Re={Re:.0f}:")
        print(f"    Smooth sphere Cd ≈ {cd_smooth_lit:.2f}  (our result: {Cd_s:.4f})")
        if cd_golf_lit:
            print(f"    Golf ball Cd   ≈ {cd_golf_lit:.2f}  (our result: {Cd_g:.4f})")
    print("="*55 + "\n")

    # ── Plot ────────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f'Golf Ball vs Smooth Sphere — Re = {Re:.0f}',
                 fontsize=15, fontweight='bold')
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    colors = {'smooth': '#2196F3', 'golf': '#FF5722'}

    # ── Panel 1: Cd time series ────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, :])

    ax1.plot(smooth['t'], smooth['Cd'],
             alpha=0.25, color=colors['smooth'], linewidth=0.8)
    ax1.plot(golf['t'], golf['Cd'],
             alpha=0.25, color=colors['golf'], linewidth=0.8)

    # Smoothed version on top
    if len(smooth['t']) > 40:
        ax1.plot(smooth['t'], moving_average(smooth['Cd']),
                 color=colors['smooth'], linewidth=2.0,
                 label=f'Smooth sphere  Cd = {Cd_s:.3f} ± {stats_s["std"]:.3f}')
        ax1.plot(golf['t'], moving_average(golf['Cd']),
                 color=colors['golf'], linewidth=2.0,
                 label=f'Golf ball      Cd = {Cd_g:.3f} ± {stats_g["std"]:.3f}')
    else:
        ax1.plot(smooth['t'], smooth['Cd'],
                 color=colors['smooth'], linewidth=2.0,
                 label=f'Smooth sphere  Cd = {Cd_s:.3f}')
        ax1.plot(golf['t'], golf['Cd'],
                 color=colors['golf'], linewidth=2.0,
                 label=f'Golf ball      Cd = {Cd_g:.3f}')

    # Mark averaging start
    ax1.axvline(t_stat, color='gray', linestyle='--', alpha=0.7,
                label=f't = {t_stat:.0f}  (averaging starts)')
    ax1.axhline(Cd_s, color=colors['smooth'], linestyle=':', alpha=0.6, linewidth=1.5)
    ax1.axhline(Cd_g, color=colors['golf'],   linestyle=':', alpha=0.6, linewidth=1.5)

    ax1.set_xlabel('Non-dimensional time  t·U₀/D', fontsize=11)
    ax1.set_ylabel('Drag coefficient  Cd', fontsize=11)
    ax1.set_title('Drag coefficient time series', fontsize=12)
    ax1.legend(fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)

    # ── Panel 2: Pressure vs viscous breakdown ─────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    categories = ['Pressure (form)\ndrag', 'Viscous (friction)\ndrag', 'Total\nCd']
    x = np.arange(len(categories))
    width = 0.32

    smooth_vals = [stats_s['Cd_p'], stats_s['Cd_v'], Cd_s]
    golf_vals   = [stats_g['Cd_p'], stats_g['Cd_v'], Cd_g]

    b1 = ax2.bar(x - width/2, smooth_vals, width,
                 label='Smooth sphere', color=colors['smooth'], alpha=0.85, edgecolor='white')
    b2 = ax2.bar(x + width/2, golf_vals,   width,
                 label='Golf ball',     color=colors['golf'],   alpha=0.85, edgecolor='white')

    ax2.bar_label(b1, fmt='%.3f', fontsize=8.5, padding=2)
    ax2.bar_label(b2, fmt='%.3f', fontsize=8.5, padding=2)

    ax2.set_ylabel('Cd component', fontsize=11)
    ax2.set_title('Drag breakdown', fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(categories, fontsize=9)
    ax2.legend(fontsize=9)
    ax2.grid(True, axis='y', alpha=0.3)
    ax2.set_ylim(bottom=0)

    # ── Panel 3: Drag reduction summary ───────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')

    summary_lines = [
        f"Re = {Re:.0f}",
        "",
        f"Smooth sphere:",
        f"  Cd = {Cd_s:.4f} ± {stats_s['std']:.4f}",
        f"  Pressure : {stats_s['Cd_p']:.4f}",
        f"  Viscous  : {stats_s['Cd_v']:.4f}",
        "",
        f"Golf ball:",
        f"  Cd = {Cd_g:.4f} ± {stats_g['std']:.4f}",
        f"  Pressure : {stats_g['Cd_p']:.4f}",
        f"  Viscous  : {stats_g['Cd_v']:.4f}",
        "",
        f"Drag change: {drag_reduction:+.1f}%",
        "",
        f"{'✓ DRAG REDUCTION' if drag_reduction > 0 else '✗ DRAG INCREASE'}",
        f"(expected above Re~4×10⁴)",
    ]

    color_result = '#2e7d32' if drag_reduction > 0 else '#c62828'
    y_pos = 0.97
    for i, line in enumerate(summary_lines):
        weight = 'bold' if (line.startswith('Re') or '±' in line or
                            'DRAG' in line or line.startswith('Golf') or
                            line.startswith('Smooth')) else 'normal'
        color = color_result if 'DRAG' in line else 'black'
        fontsize = 11 if 'DRAG' in line else 9.5
        ax3.text(0.05, y_pos - i*0.066, line,
                 transform=ax3.transAxes,
                 fontsize=fontsize, fontweight=weight,
                 color=color, fontfamily='monospace')

    ax3.set_title('Summary', fontsize=12)
    rect = plt.Rectangle((0, 0), 1, 1, fill=False,
                          edgecolor='gray', linewidth=1.5,
                          transform=ax3.transAxes)
    ax3.add_patch(rect)

    # ── Save ────────────────────────────────────────────────────────────────────
    outfile = 'drag_comparison.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"Plot saved: {outfile}")

    # Also write a machine-readable CSV summary
    with open('drag_summary.csv', 'w') as f:
        f.write('case,Re,Cd_mean,Cd_std,Cd_pressure,Cd_viscous\n')
        f.write(f"smooth,{Re},{Cd_s:.6f},{stats_s['std']:.6f},"
                f"{stats_s['Cd_p']:.6f},{stats_s['Cd_v']:.6f}\n")
        f.write(f"golf,{Re},{Cd_g:.6f},{stats_g['std']:.6f},"
                f"{stats_g['Cd_p']:.6f},{stats_g['Cd_v']:.6f}\n")
        f.write(f"drag_reduction_pct,{Re},{drag_reduction:.4f},,,\n")
    print("Summary CSV saved: drag_summary.csv")


if __name__ == "__main__":
    main()
