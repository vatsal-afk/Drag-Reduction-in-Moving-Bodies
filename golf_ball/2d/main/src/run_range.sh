#!/bin/bash
# run_range.sh - Run Basilisk simulation over a range of Reynolds numbers
set -e

# Compile the code
echo "Compiling..."
qcc -O2 -Wall -DBUMPY test_final.c -o test_final -lm

# Array of Reynolds numbers (laminar to turbulent)
RE_LIST=("100" "500" "1000" "5000" "50000" "100000")

# Organized structure in results/
BASE_DIR="../results/parametric_sweep"
mkdir -p "$BASE_DIR"

echo "Starting parameter sweep..."
for RE in "${RE_LIST[@]}"; do
    echo "==========================================="
    echo " Running simulation for Re = $RE "
    echo "==========================================="
    
    # Run the simulation (add mpirun if you want parallel)
    ./test_final "$RE"
    
    # Move results to their own folder so vtu files don't clutter everything
    DIR="$BASE_DIR/Re_$RE"
    mkdir -p "$DIR"
    mv drag_*_Re"$RE".dat "$DIR/" 2>/dev/null || true
    mv vort_*_Re"$RE".mp4 "$DIR/" 2>/dev/null || true
    mv *_Re"$RE"-*.vtu "$DIR/" 2>/dev/null || true
done

echo "Done! Results are separated into $BASE_DIR"
echo "Run: python3 plot_sweep.py to view the analysis"

