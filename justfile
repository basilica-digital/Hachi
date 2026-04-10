
run:
    RAYON_NUM_THREADS=16 RUSTFLAGS="-C target-cpu=native -C codegen-units=1" cargo run --release
