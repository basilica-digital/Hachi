use crate::hachi::prove::SumcheckProof;
use crate::sumcheck::verify::sumcheck_verify;
use crate::hachi::setup::SetupParams;
use crate::utils::random::{
    generate_random_eval_points_q_seeded, generate_random_q_element_seeded,
    generate_sparse_c_idx_seeded,
};
use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;

#[derive(Debug, Clone)]
pub struct Challenge{
    pub a0: Vec<u32>,
    pub a1: Vec<u32>,
    pub a2: Vec<u32>,
    pub a3: Vec<u32>,
    pub b0: Vec<u32>,
    pub b1: Vec<u32>,
    pub b2: Vec<u32>,
    pub b3: Vec<u32>,
    pub c: Vec<i16>,
    pub alpha: Vec<u32>,
    pub tau_0: [[u32; 4]; 28],
    pub tau_1: [[u32; 4]; 4]
}

/// Deterministically derive the interactive challenge from `seed`.
///
/// NOTE: this is a stop-gap so the split `prove` / `verify` CLI binaries
/// produce matching challenges. It is **not** cryptographically sound —
/// the prover can pick the seed. Replace with Fiat-Shamir (hash of
/// commitment + public params) before relying on the proof's soundness.
pub fn sample_challenge(params: &SetupParams, seed: u64) -> Challenge {
    let mut rng = ChaCha12Rng::seed_from_u64(seed);
    let a0 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let a1 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let a2 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let a3 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let b0 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let b1 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let b2 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let b3 = generate_random_eval_points_q_seeded(params.n, 1, &mut rng);
    let c = generate_sparse_c_idx_seeded(params.n, &mut rng);
    let alpha = generate_random_q_element_seeded(1, 16, &mut rng);
    let mut tau_0 = [[0u32; 4]; 28];
    let mut tau_1 = [[0u32; 4]; 4];
    for i in 0..28 {
        if i < 4 {
            tau_1[i] = generate_random_q_element_seeded(1, 16, &mut rng)[0..4]
                .try_into()
                .unwrap();
        }
        tau_0[i] = generate_random_q_element_seeded(1, 16, &mut rng)[0..4]
            .try_into()
            .unwrap();
    }
    Challenge {
        a0: a0.to_vec(),
        a1: a1.to_vec(),
        a2: a2.to_vec(),
        a3: a3.to_vec(),
        b0: b0.to_vec(),
        b1: b1.to_vec(),
        b2: b2.to_vec(),
        b3: b3.to_vec(),
        c: c.to_vec(),
        alpha: alpha.to_vec(),
        tau_0,
        tau_1,
    }
}

pub fn verify(
    params: &SetupParams,
    proof: &SumcheckProof,
    target_sum_alg: [u32; 4],
    challenge: &Challenge
) -> bool {
    let is_valid = unsafe {sumcheck_verify(proof, params, target_sum_alg, [0u32, 0, 0, 0], &params.w_i, challenge)};
    is_valid
}

