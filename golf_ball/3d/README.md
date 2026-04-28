# Golf Ball CFD Simulation (Basilisk)

This repository contains a 3D Computational Fluid Dynamics (CFD) simulation comparing the aerodynamic drag of a dimpled golf ball against a smooth sphere at a specific Reynolds number. 

The simulation is built on the [Basilisk](http://basilisk.fr/) framework, utilizing an adaptive octree mesh and embedded boundaries to resolve complex micro-geometry (dimples) efficiently.

## 🚀 Setup Instructions for a New System

Since Basilisk uses a specialized C compiler (`qcc`), there are a few system-level dependencies required before you can run the pipeline.

### 1. Install Basilisk

Basilisk is compiled from source. Run these commands on your new machine:

```bash
# Get required system dependencies
sudo apt-get update
sudo apt-get install gcc make gawk gnuplot

# Download and compile Basilisk
cd ~
wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar -xzf basilisk.tar.gz
cd basilisk/src
make

# Add qcc to your PATH (You should add these lines to your ~/.bashrc)
export BASILISK=$HOME/basilisk/src
export PATH=$PATH:$BASILISK
```

### 2. Set up the Python Environment

The Makefile relies on a local Python virtual environment to cleanly repair the 3D STL geometry and generate the final drag comparison graphs.

```bash
# Inside the golf-ball project directory:
python3 -m venv venv
./venv/bin/pip install -r requirements.txt
```

### 3. Running the Pipeline

Once Basilisk and Python are set up, you can run the entire pipeline end-to-end with a single command:

```bash
make all
```

This will automatically:
1. Repair and normalize the raw STL file to a clean binary format (`step1_preprocess_stl.py`).
2. Compile and run the `smooth_sphere` simulation.
3. Compile and run the `golf_ball` simulation using a highly optimized pre-refinement strategy.
4. Process the drag data and output comparison charts (`step4_analyze.py`).

## ⚙️ Advanced Configuration (MPI Parallelization)

By default, the simulation runs on a single CPU core. 3D CFD is highly computationally expensive. If your system has multiple cores, you can significantly speed up the simulation using OpenMPI.

1. Install MPI on your system:
   ```bash
   sudo apt-get install openmpi-bin libopenmpi-dev
   ```
2. Open the `Makefile` and locate the **Compiler flags** and **Configuration** sections.
3. Change `MPI_RANKS := 1` to the number of CPU cores you wish to dedicate (e.g., `MPI_RANKS := 8`).
4. Comment out the standard `QCC_FLAGS` and uncomment the MPI-enabled ones:
   ```makefile
   # QCC_FLAGS := -O2 -Wall -disable-dimensions
   # CC        := mpicc -D_MPI=$(MPI_RANKS)
   ```

## 🎥 Visualization (ParaView)

During execution, both simulations will output `.vtu` (VTK Unstructured Grid) snapshots every 5 simulation seconds into `vtu_smooth/` and `vtu_golf/` respectively. 

You can load these directories into [ParaView](https://www.paraview.org/) to create 3D renders and videos of the wake turbulence (vorticity/pressure) and the geometric surface (`cs` fraction = 0.5).

*Note: These directories are ignored by Git due to their massive size.*
