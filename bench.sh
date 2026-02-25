#!/bin/bash
# THREAD_COUNTS=(1 2 4 8 16 32)
THREAD_COUNTS=(16)
EXECUTABLE="./target/release/ntt_avx512_test"
echo "========================================"
for threads in "${THREAD_COUNTS[@]}"; do
    echo "▶ RAYON_NUM_THREADS=$threads"
    RAYON_NUM_THREADS=$threads $EXECUTABLE
    echo "========================================"
done