//! Full commit / prove / verify pipeline coverage.
//!
//! These are the regression tests that actually exercise the SIMD cores:
//! the happy path returns a proof that verifies, and a handful of tamper
//! variants ensure `verify` doesn't silently accept garbage. They are
//! marked `#[ignore]` because each run allocates ~4 GiB for the witness
//! plus several GiB of intermediate prover state and takes minutes.
//!
//! Run them with:
//!
//! ```sh
//! RUSTFLAGS="-C target-cpu=native" cargo test --release --test end_to_end -- --ignored --test-threads=1
//! ```
//!
//! or via `just test`. `--test-threads=1` matters: the fixture is ~10 GiB
//! at peak and running tests in parallel OOMs quickly.
//!
//! AVX-512 is required; the crate's intrinsics assume it.
#![allow(non_snake_case)]

use std::sync::OnceLock;

use hachi::hachi::commit::Commit;
use hachi::hachi::prove::{prove, SumcheckProof};
use hachi::hachi::setup::{setup_with_seed, SetupParams};
use hachi::hachi::verify::{sample_challenge, verify, Challenge};
use hachi::utils::ds::AlignedU8Vec;
use hachi::utils::random::generate_random_bytes_4bit_packed;

const SEED: u64 = 0;

/// One-shot pipeline output shared across all tamper tests. Building it
/// costs minutes and several GiB, so we pay that price exactly once per
/// `cargo test` invocation even when multiple `#[test]`s touch it.
struct Artifacts {
    params: SetupParams,
    challenge: Challenge,
    proof: SumcheckProof,
    target_sum: [u32; 4],
}

fn artifacts() -> &'static Artifacts {
    static CELL: OnceLock<Artifacts> = OnceLock::new();
    CELL.get_or_init(build_artifacts)
}

fn build_artifacts() -> Artifacts {
    let params = setup_with_seed(SEED);

    // Witness bytes are drawn from the OS RNG — not seeded — so this
    // test also catches bugs that only show up for generic random
    // witnesses rather than a specific fixed pattern.
    let witness_bytes: usize = params
        .witness_size_bytes()
        .try_into()
        .expect("witness size exceeds usize");
    let witness: AlignedU8Vec = generate_random_bytes_4bit_packed(witness_bytes);

    let challenge = sample_challenge(&params, SEED);
    let commitment = Commit(&params, &witness);
    let (proof, target_sum) =
        unsafe { prove(&params, &witness, commitment, &challenge) };

    Artifacts {
        params,
        challenge,
        proof,
        target_sum,
    }
}

#[test]
#[ignore = "full pipeline: ~minutes runtime, ~10 GiB peak RAM, needs AVX-512"]
fn commit_prove_verify_happy_path() {
    let a = artifacts();
    assert!(
        verify(&a.params, &a.proof, a.target_sum, &a.challenge),
        "honest proof must verify"
    );
}

#[test]
#[ignore = "full pipeline: ~minutes runtime, ~10 GiB peak RAM, needs AVX-512"]
fn verify_rejects_tampered_target_sum() {
    let a = artifacts();
    let mut bad = a.target_sum;
    bad[0] = bad[0].wrapping_add(1);
    assert!(
        !verify(&a.params, &a.proof, bad, &a.challenge),
        "verify must reject a corrupted target sum"
    );
}

#[test]
#[ignore = "full pipeline: ~minutes runtime, ~10 GiB peak RAM, needs AVX-512"]
fn verify_rejects_tampered_round_polynomial() {
    let a = artifacts();
    // Bumping the first coefficient of the first X-round polynomial
    // breaks the `p_alg[0] + p_alg[1] == target_sum` check immediately.
    let mut bad = a.proof.clone();
    bad.p_alg_x[0][0][0] = bad.p_alg_x[0][0][0].wrapping_add(1);
    assert!(
        !verify(&a.params, &bad, a.target_sum, &a.challenge),
        "verify must reject a tampered round polynomial"
    );
}

#[test]
#[ignore = "full pipeline: ~minutes runtime, ~10 GiB peak RAM, needs AVX-512"]
fn verify_rejects_tampered_final_w() {
    let a = artifacts();
    let mut bad = a.proof.clone();
    bad.final_w[0] = bad.final_w[0].wrapping_add(1);
    assert!(
        !verify(&a.params, &bad, a.target_sum, &a.challenge),
        "verify must reject a tampered final evaluation"
    );
}

#[test]
#[ignore = "full pipeline: ~minutes runtime, ~10 GiB peak RAM, needs AVX-512"]
fn verify_rejects_wrong_challenge() {
    let a = artifacts();
    // Same params, different challenge seed: all of tau/alpha/c flip, so
    // the sumcheck's first round target no longer matches the claimed
    // sum the prover baked in.
    let wrong = sample_challenge(&a.params, SEED.wrapping_add(1));
    assert!(
        !verify(&a.params, &a.proof, a.target_sum, &wrong),
        "verify must reject a proof bound to a different challenge"
    );
}
