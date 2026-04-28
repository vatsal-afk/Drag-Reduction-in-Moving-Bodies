# Makefile — Golf Ball CFD Drag Comparison
# ==========================================
# Prerequisites:
#   - Basilisk installed (http://basilisk.fr/src/INSTALL)
#   - qcc available in PATH
#   - Python 3 with trimesh, numpy, matplotlib
#
# Usage:
#   make preprocess   — clean and normalize the STL (step 1)
#   make smooth       — compile and run smooth sphere (step 2)
#   make golf         — compile and run golf ball (step 3)
#   make analyze      — compare results and plot (step 4)
#   make all          — run full pipeline
#   make clean        — remove binaries and data

# ── Configuration ──────────────────────────────────────────────────────────────

# Reynolds number (change here to run at different Re)
RE        := 50000.0

# Refinement levels
GEOM_LEVEL := 10
FLOW_LEVEL := 10

# Number of MPI ranks (set to 1 for serial)
MPI_RANKS := 1

# Raw STL from Sketchfab download
RAW_STL   := golf_ball_raw.stl

# Preprocessed STL (input to golf_ball.c)
GOLF_STL  := golfball_clean.stl

# ── Compiler flags ──────────────────────────────────────────────────────────────

# -disable-dimensions: suppresses Basilisk's compile-time dimensional analysis
# (required by viscosity-embed.h — separate from the GL issue which is now fixed)
QCC_FLAGS := -O2 -Wall -disable-dimensions

# For MPI parallel, uncomment this and comment out the above:
# QCC_FLAGS := -O2 -Wall -disable-dimensions
# CC        := mpicc -D_MPI=$(MPI_RANKS)

# ── Targets ────────────────────────────────────────────────────────────────────

.PHONY: all preprocess smooth golf analyze clean help

all: preprocess smooth golf analyze

help:
	@echo ""
	@echo "Golf Ball CFD Drag Comparison"
	@echo "=============================="
	@echo ""
	@echo "Targets:"
	@echo "  make preprocess    Clean/normalize STL (needs $(RAW_STL))"
	@echo "  make smooth        Compile + run smooth sphere simulation"
	@echo "  make golf          Compile + run golf ball simulation"
	@echo "  make analyze       Plot drag comparison"
	@echo "  make all           Full pipeline"
	@echo ""
	@echo "Configuration (edit Makefile or override):"
	@echo "  RE=$(RE)           Reynolds number"
	@echo "  GEOM_LEVEL=$(GEOM_LEVEL)     Geometry refinement level"
	@echo "  FLOW_LEVEL=$(FLOW_LEVEL)      Flow AMR level"
	@echo "  MPI_RANKS=$(MPI_RANKS)        MPI processes (1=serial)"
	@echo ""

# ── Step 1: Preprocess STL ─────────────────────────────────────────────────────

preprocess: $(GOLF_STL)

$(GOLF_STL): $(RAW_STL) step1_preprocess_stl.py
	@echo ""
	@echo "=== Step 1: Preprocessing STL ==="
	./venv/bin/python3 step1_preprocess_stl.py $(RAW_STL) $(GOLF_STL)

# ── Step 2: Smooth Sphere ──────────────────────────────────────────────────────

smooth: smooth_sphere drag_smooth.dat

smooth_sphere: smooth_sphere.c embed_force_3d.h Makefile
	@echo ""
	@echo "=== Step 2: Compiling smooth_sphere ==="
	qcc $(QCC_FLAGS) \
	    -DRE=$(RE) \
	    -DMAXLEVEL=$(FLOW_LEVEL) \
	    smooth_sphere.c -o smooth_sphere -lm
	@echo "Compiled: smooth_sphere"

drag_smooth.dat: smooth_sphere
	@echo ""
	@echo "=== Step 2: Running smooth sphere (Re=$(RE)) ==="
	@echo "This may take 10-60 minutes depending on Re and level."
	@mkdir -p vtu_smooth
ifeq ($(MPI_RANKS),1)
	./smooth_sphere 2>&1 | tee smooth_sphere.log
else
	mpirun -n $(MPI_RANKS) ./smooth_sphere 2>&1 | tee smooth_sphere.log
endif
	@echo "Done. Output: drag_smooth.dat"

# ── Step 3: Golf Ball ──────────────────────────────────────────────────────────

golf: golf_ball drag_golf.dat

golf_ball: golf_ball.c embed_force_3d.h $(GOLF_STL) Makefile
	@echo ""
	@echo "=== Step 3: Compiling golf_ball ==="
	qcc $(QCC_FLAGS) \
	    -DRE=$(RE) \
	    -DGEOM_LEVEL=$(GEOM_LEVEL) \
	    -DFLOW_LEVEL=$(FLOW_LEVEL) \
	    -DSTL_FILE=$(GOLF_STL) \
	    golf_ball.c -o golf_ball -lm
	@echo "Compiled: golf_ball"

drag_golf.dat: golf_ball $(GOLF_STL)
	@echo ""
	@echo "=== Step 3: Running golf ball (Re=$(RE)) ==="
	@echo "Geometry refinement to level $(GEOM_LEVEL) may take extra time."
	@mkdir -p vtu_golf
ifeq ($(MPI_RANKS),1)
	./golf_ball 2>&1 | tee golf_ball.log
else
	mpirun -n $(MPI_RANKS) ./golf_ball 2>&1 | tee golf_ball.log
endif
	@echo "Done. Output: drag_golf.dat"

# ── Step 4: Analysis ──────────────────────────────────────────────────────────

analyze: drag_comparison.png

drag_comparison.png: drag_smooth.dat drag_golf.dat step4_analyze.py
	@echo ""
	@echo "=== Step 4: Analyzing results ==="
	./venv/bin/python3 step4_analyze.py
	@echo "Done. Output: drag_comparison.png, drag_summary.csv"

# ── Clean ─────────────────────────────────────────────────────────────────────

clean:
	rm -f smooth_sphere golf_ball
	rm -f drag_smooth.dat drag_golf.dat
	rm -f result_smooth.txt result_golf.txt
	rm -f drag_comparison.png drag_summary.csv
	rm -f *.log *.mp4 *.ppm
	rm -rf vtu_smooth vtu_golf
	@echo "Cleaned."

cleanall: clean
	rm -f $(GOLF_STL)