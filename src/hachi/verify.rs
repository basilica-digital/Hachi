use crate::hachi::prove::SumcheckProof;
use crate::sumcheck::verify::sumcheck_verify;
use crate::hachi::setup::SetupParams;
use crate::utils::random::generate_random_eval_points_q;
use crate::utils::random::generate_sparse_c_idx;
use crate::utils::random::generate_random_q_element;

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

pub fn sample_challenge(params: &SetupParams) -> Challenge {
    let a0 = generate_random_eval_points_q(params.n, 1);
    let a1 = generate_random_eval_points_q(params.n, 1);
    let a2 = generate_random_eval_points_q(params.n, 1);
    let a3 = generate_random_eval_points_q(params.n, 1);
    let b0 = generate_random_eval_points_q(params.n, 1);
    let b1 = generate_random_eval_points_q(params.n, 1);
    let b2 = generate_random_eval_points_q(params.n, 1);
    let b3 = generate_random_eval_points_q(params.n, 1);
    let c = generate_sparse_c_idx(params.n);
    let alpha = generate_random_q_element(1, 16);
    let mut tau_0 = [[0u32; 4]; 28];
    let mut tau_1 = [[0u32; 4]; 4];
    for i in 0..28{
        if i<4{
            tau_1[i] = generate_random_q_element(1, 16)[0..4].try_into().unwrap();
        }
        tau_0[i] = generate_random_q_element(1, 16)[0..4].try_into().unwrap();
    }
    Challenge { a0:a0.to_vec(), a1:a1.to_vec(), a2:a2.to_vec(), a3:a3.to_vec(), b0:b0.to_vec(), b1:b1.to_vec(), b2:b2.to_vec(), b3:b3.to_vec(), c:c.to_vec(), alpha:alpha.to_vec(), tau_0, tau_1 }
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

