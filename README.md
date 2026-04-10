# Hachi
A lattice-based polynomial commitment scheme

## Overview

Hachi implements a sumcheck-based proof system over a lattice commitment of the form `u = As - Bt`. The scheme operates over the extension field `GF((2^32 - 99)^4)` defined by the irreducible polynomial `x^4 - 2`.

The prover runs a 27-round sumcheck (17 "X" rounds + 10 "Y" rounds) with two constraint types:
- **Algebraic constraints** — degree-2 polynomials (3 evaluation points per round)
- **Norm/balance constraints** — degree-18 polynomials (19 evaluation points per round) encoding bounded-norm checks via `c_bal = (w^2 + 8w) * prod(w^2 - k^2)` for `k = 1..7`

## Parameters

| Parameter | Value |
|-----------|-------|
| Ring dimension (`n`) | 2^10 = 1024 |
| Matrix height | 2^13 = 8192 |
| Constraints | 11 |
| Field modulus (`q`) | 2^32 - 99 |
| NTT primes (RNS) | 759207937, 759304193, 2079301633, 2079305729 |
| Fiat-Shamir hash | BLAKE3 -> ChaCha12 |

## Dependencies

- Rust (edition 2021)
- x86-64 CPU with AVX-512 support
- [just](https://github.com/casey/just) (command runner)

## Building and running

```
just run
```

This compiles in release mode with `target-cpu=native` and runs with 16 Rayon threads. Equivalent to:

```
RAYON_NUM_THREADS=16 RUSTFLAGS="-C target-cpu=native -C codegen-units=1" cargo run --release
```

## Project structure

```
src/
├── main.rs              # Entry point: setup -> commit -> prove -> verify
├── math/
│   ├── fields.rs        # Scalar GF(p^4) arithmetic
│   ├── field_simd.rs    # AVX-512 vectorized field ops (Barrett, Montgomery, 8-wide extension mul)
│   ├── rns.rs           # RNS decomposition/reconstruction for NTT primes
│   └── eq.rs            # Equality polynomial table construction
├── hachi/
│   ├── setup.rs         # Public parameter generation
│   ├── commit.rs        # Lattice commitment (u = As - Bt)
│   ├── prove.rs         # Prover: witness preparation + sumcheck invocation
│   └── verify.rs        # Verifier
├── sumcheck/
│   ├── prove.rs         # 27-round sumcheck prover
│   └── verify.rs        # Sumcheck verifier
├── mlp/
│   ├── mle.rs           # Multilinear extensions for witness W and constraint matrix M
│   └── eval.rs          # Verifier-side evaluation (eq, alpha, step indicators)
└── utils/
    ├── ds.rs            # 64-byte aligned vector types for SIMD
    ├── random.rs        # Random data generation and Fiat-Shamir challenges
    └── fs.rs            # BLAKE3-based Fiat-Shamir transcript
```

## License

Apache-2.0
