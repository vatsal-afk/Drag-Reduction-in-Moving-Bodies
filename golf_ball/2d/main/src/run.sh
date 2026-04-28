#!/bin/bash

# 1. Clean up old data
rm -f drag_smooth.dat drag_bumpy.dat *.vtu.vtu

# 2. Compile Smooth version
qcc -O2 cyl2d.c -o smooth_sim -lm

# 3. Compile Bumpy version (with -DBUMPY flag)
qcc -O2 -DBUMPY cyl2d.c -o bumpy_sim -lm

echo "Starting simulations... this may take some time."

# 4. Run both in parallel (if you have MPI)
# If you DON'T have MPI, remove "-D_MPI=1" from compilation 
# and just run "./smooth_sim & ./bumpy_sim &"
./smooth_sim > log_smooth 2>&1 &
./bumpy_sim > log_bumpy 2>&1 &

echo "Both simulations are running in the background."
echo "Use 'tail -f log_smooth' to watch progress."

wait
echo "Simulations complete! Use your plot.py to see the results."