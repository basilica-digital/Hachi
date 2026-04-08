use crate::hachi::prove::SumcheckProof;
use crate::sumcheck::verify::sumcheck_verify;
use crate::hachi::setup::SetupParams;

pub fn verify(
    params: &SetupParams,
    proof: &SumcheckProof,
    target_sum_alg: [u32; 4],
    tau_0: &[[u32; 4]; 27],
    alpha: &[u32; 4],
) -> bool {
    let is_valid = sumcheck_verify(proof, target_sum_alg, [0u32, 0, 0, 0], &params.w_i, tau_0, alpha);
    is_valid
}