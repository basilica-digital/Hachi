# Hachi
A lattice-based polynomial commitment scheme

Test commit (``src/commit.rs``):
```
RUSTFLAGS="-C target-cpu=native -C codegen-units=1" cargo build --release
RAYON_NUM_THREADS=1 ./target/release/Hachi30
```

