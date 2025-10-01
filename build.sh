#!/bin/bash

# build.sh: Main control script for the signal analysis project.

# --- Configuration ---
CSV_FILE="mental-state.csv"
TARGET_COLUMN="lag1_mean_0"
EXECUTABLE="./build/analyzer"
PLOTTER_SCRIPT="scripts/plotter.py"

# Output directories
RESULTS_DIR="results"
DATA_DIR="$RESULTS_DIR/data"
PLOTS_DIR="$RESULTS_DIR/plots"
PERF_DIR="$RESULTS_DIR/performance"

# --- Functions ---

clean_up() {
    echo "--- Cleaning Up Previous Build and Results ---"
    make clean > /dev/null 2>&1
    rm -rf "$RESULTS_DIR"
    echo "Cleanup complete."
    echo ""
}

command_exists() { command -v "$1" >/dev/null 2>&1; }

check_prerequisites() {
    echo "--- Checking Prerequisites ---"
    if ! command_exists gcc || ! command_exists make; then
        echo "Error: GCC and Make are required." >&2; exit 1;
    fi
    if ! command_exists python3; then
        echo "Error: Python 3 is required." >&2; exit 1;
    fi
    # Added scipy to the check
    python3 -c "import numpy, matplotlib, scipy" 2>/dev/null
    if [ $? -ne 0 ]; then
        echo "Error: Python libraries numpy, matplotlib, and scipy are required ('pip install numpy matplotlib scipy')." >&2; exit 1;
    fi
    if [ ! -f /usr/include/fftw3.h ] && [ ! -f /usr/local/include/fftw3.h ]; then
        echo "Warning: FFTW3 header not found. Please install it (e.g., 'sudo apt-get install libfftw3-dev')."
    fi
    echo "Prerequisites met."
    echo ""
}

build_c_code() {
    echo "--- Compiling C Code ---"
    make
    if [ $? -ne 0 ]; then
        echo "Error: Compilation failed. Aborting."
        exit 1
    fi
    echo ""
}

run_analysis() {
    echo "--- Running C Analyzer ---"
    if [ ! -f "$CSV_FILE" ]; then
        echo "Error: Input file '$CSV_FILE' not found."
        exit 1
    fi
    
    mkdir -p "$DATA_DIR" "$PLOTS_DIR" "$PERF_DIR"
    
    $EXECUTABLE "$CSV_FILE" "$TARGET_COLUMN" "$DATA_DIR"
    if [ $? -ne 0 ]; then
        echo "Error: C analyzer failed. Aborting."
        exit 1
    fi
    echo "C analysis complete. Data files generated in '$DATA_DIR'."
    echo ""
}

generate_plots() {
    echo "--- Generating Plots ---"
    echo "Plots will be shown in pop-up windows. Close all windows to proceed to the next."
    python3 $PLOTTER_SCRIPT
    if [ $? -ne 0 ]; then
        echo "Error: Python plotting script failed."
        exit 1
    fi
    echo "Plots also saved in '$PLOTS_DIR'."
    echo ""
}

# --- Main Logic ---
case "$1" in
    clean)
        clean_up
        ;;
    *)
        clean_up
        check_prerequisites
        build_c_code
        run_analysis
        generate_plots
        echo "--- Project Execution Finished Successfully ---"
        ;;
esac

exit 0