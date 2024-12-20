#!/bin/bash

export LC_ALL=C

# Usage check
if [ $# -lt 1 ]; then
    echo "Usage: $0 results_directory"
    exit 1
fi

RESULTS_DIR="$1"

# Create results directory if it doesn't exist
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Creating directory $RESULTS_DIR"
    mkdir -p "$RESULTS_DIR"
fi

# Check if directory creation was successful
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Directory $RESULTS_DIR does not exist and could not be created"
    exit 2
fi

# Find the highest existing numeric directory name
N=$( ( echo 0; ls -1 "$RESULTS_DIR" ) | grep -E '^[0-9]+$' | sort -n | tail -n 1 | awk '{print $1}')

NN=0
while [ "$NN" == 0 ]; do
    N=$(( N + 1 ))
    if [ $N -ge 1000 ]; then
        echo "Too busy directory"
        exit 3
    fi
    N3=$( printf "%03d" "$N" )
    if [ ! -e "$RESULTS_DIR/$N3" ] && [ ! -e "$RESULTS_DIR/${N3}_big2" ]; then
        NN="$N3"
    fi
done

echo "Starting from $NN"

# Check for Intel MPI compiler
OPT=""
if which mpif90 | grep -q "intel"; then
    OPT="-c benchmarking_ifx"
fi

# Run benchmarks
for i in $( seq $N $(( N + 2 )) ); do
    ./benchmarking/piernik_bench.sh ${OPT} | tee "$RESULTS_DIR/$( printf "%03d" "$i" )"
done

for b in 1.5 2; do
    BIG=$b ./benchmarking/piernik_bench.sh ${OPT} | tee "$RESULTS_DIR/$( printf "%03d" "$N" )_big$b"
done
