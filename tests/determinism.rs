//! Fast regression tests that don't allocate the full witness.
//!
//! These cover the deterministic pieces of the pipeline — public
//! parameters and the Fiat-Shamir-ish challenge derivation — so the
//! default `cargo test` run catches accidental reshuffles of the seeded
//! RNG streams without paying the cost of a full `commit`/`prove`.
//!
//! Requires an AVX-512 host (the crate's intrinsics assume it). Build
//! with `RUSTFLAGS="-C target-cpu=native"` — see `justfile`'s
//! `test-fast` recipe.
//!
//! Slow end-to-end coverage lives in `tests/end_to_end.rs`.
#![allow(non_snake_case)]

use hachi::hachi::setup::{setup_with_seed, SetupParams};
use hachi::hachi::verify::{sample_challenge, Challenge};

const SEED: u64 = 0;

/// Assert that every field we can cheaply compare matches between two
/// setups. `plans` is skipped — `tfhe_ntt::Plan` isn't `PartialEq` and
/// it's a pure function of the moduli, which are hard-coded.
fn assert_setup_eq(a: &SetupParams, b: &SetupParams) {
    assert_eq!(a.constraints, b.constraints, "constraints");
    assert_eq!(a.height_2, b.height_2, "height_2");
    assert_eq!(a.height_4, b.height_4, "height_4");
    assert_eq!(a.n, b.n, "n");
    assert_eq!(a.q, b.q, "q");
    assert_eq!(a.q64, b.q64, "q64");
    assert_eq!(a.g_matrix_2, b.g_matrix_2, "g_matrix_2");
    assert_eq!(a.g_matrix_4, b.g_matrix_4, "g_matrix_4");
    assert_eq!(a.w_i, b.w_i, "w_i");
    assert_eq!(&*a.d, &*b.d, "d");
    assert_eq!(&*a.e, &*b.e, "e");
}

fn assert_challenge_eq(a: &Challenge, b: &Challenge) {
    assert_eq!(a.a0, b.a0, "a0");
    assert_eq!(a.a1, b.a1, "a1");
    assert_eq!(a.a2, b.a2, "a2");
    assert_eq!(a.a3, b.a3, "a3");
    assert_eq!(a.b0, b.b0, "b0");
    assert_eq!(a.b1, b.b1, "b1");
    assert_eq!(a.b2, b.b2, "b2");
    assert_eq!(a.b3, b.b3, "b3");
    assert_eq!(a.c, b.c, "c");
    assert_eq!(a.alpha, b.alpha, "alpha");
    assert_eq!(a.tau_0, b.tau_0, "tau_0");
    assert_eq!(a.tau_1, b.tau_1, "tau_1");
}

#[test]
fn setup_is_deterministic_for_same_seed() {
    let a = setup_with_seed(SEED);
    let b = setup_with_seed(SEED);
    assert_setup_eq(&a, &b);
}

#[test]
fn setup_differs_across_seeds() {
    let a = setup_with_seed(0);
    let b = setup_with_seed(1);
    // The seeded RNG feeds `d` and `e`; every other field is fixed by the
    // scheme parameters. Checking `d` is sufficient to confirm the seed
    // actually threads through to `generate_random_q_element_seeded`.
    assert_ne!(&*a.d, &*b.d, "d should diverge across seeds");
}

#[test]
fn sample_challenge_is_deterministic_for_same_seed() {
    let params = setup_with_seed(SEED);
    let a = sample_challenge(&params, SEED);
    let b = sample_challenge(&params, SEED);
    assert_challenge_eq(&a, &b);
}

#[test]
fn sample_challenge_differs_across_seeds() {
    let params = setup_with_seed(SEED);
    let a = sample_challenge(&params, 0);
    let b = sample_challenge(&params, 1);
    // Any divergent field is enough to prove the seed isn't being
    // ignored; pick the smallest one so failures are easy to read.
    assert_ne!(a.alpha, b.alpha, "alpha should diverge across seeds");
}

#[test]
fn setup_matches_scheme_parameters() {
    // Pins the public parameters so silently widening a matrix dimension
    // (which would invalidate existing commitments) trips a test instead
    // of a runtime shape mismatch deep inside `prove`.
    let params = setup_with_seed(SEED);
    assert_eq!(params.constraints, 15);
    assert_eq!(params.height_2, 1 << 14);
    assert_eq!(params.height_4, 1 << 13);
    assert_eq!(params.n, 1 << 10);
    assert_eq!(params.q, 4294967197);
    assert_eq!(params.q64, 4294967197);
    assert_eq!(params.witness_size_bytes(), (1u64 << 13) * (1u64 << 10) * (1u64 << 10) / 2);
}
