#!/bin/bash

# Assumes all requirements satisfied in user space

# Presets
BUILD_DIR="build"
BUILD_TARGET="Benchmark"

# Build project
cd $(pwd)/.. || exit 1

mkdir -p $BUILD_DIR
cd $BUILD_DIR || exit 1

cmake .. \
-DCMAKE_BUILD_TYPE=Release \
-DFFT_ENABLE_TESTS=OFF \
-DFFT_ENABLE_VECT_REPORT=OFF \
-DFFT_ENABLE_MANUAL_VECT=OFF

ninja $BUILD_TARGET -j0

# After successful build
if [ ! -f "projects/$BUILD_TARGET/$BUILD_TARGET" ]; then
    echo "projects/$BUILD_TARGET/$BUILD_TARGET doesn't exist" && exit 1
fi

# Set the number of nodes
SBATCH --nodes=1
# Hyper-threading off
SBATCH --threads-per-core=1
# Set max wallclock time
SBATCH --time=6-0
# Set name of job
SBATCH --job-name=fft_makarov_calculation
# Set queue name
SBATCH -p broadwell
SBATCH --ntasks-per-node=1

./$BUILD_TARGET 2>&1 | tee output.txt

# Get results

cat output.txt
