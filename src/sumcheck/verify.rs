use crate::math::fields::*;
use crate::mlp::eval::*;
use crate::Q;
use crate::hachi::prove::SumcheckProof;

pub fn sumcheck_verify(
    proof: &SumcheckProof,
    mut target_sum_alg: [u32; 4],
    mut target_sum_norm: [u32; 4],
    w_i: &[u32; 19],
    tau_0: &[[u32; 4]; 27],
    alpha: &[u32; 4]
) -> bool {
    let q = 4294967197u32;
    // X
    for round in 0..17 {
        let p_alg = &proof.p_alg_x[round];
        let p_norm = &proof.p_norm_x[round];
        let r = &proof.r_x[round]; 

        let current_sum_alg = unsafe { extadd(&p_alg[0], &p_alg[1]) };
        if current_sum_alg != target_sum_alg { return false; }
        
        let current_sum_norm = unsafe { extadd(&p_norm[0], &p_norm[1]) };
        if current_sum_norm != target_sum_norm { return false; }
        
        target_sum_alg = unsafe { eval_quadratic(&p_alg[0], &p_alg[1], &p_alg[2], r) };
        target_sum_norm = unsafe { eval_degree_18(p_norm, r, w_i) };
    }
    // Y
    for round in 0..10 {
        let p_alg = &proof.p_alg_y[round];
        let p_norm = &proof.p_norm_y[round];
        let r = &proof.r_y[round];

        let current_sum_alg = unsafe { extadd(&p_alg[0], &p_alg[1]) };
        if current_sum_alg != target_sum_alg { return false; }
        
        let current_sum_norm = unsafe { extadd(&p_norm[0], &p_norm[1]) };
        if current_sum_norm != target_sum_norm { return false; }
        
        target_sum_alg = unsafe { eval_quadratic(&p_alg[0], &p_alg[1], &p_alg[2], r) };
        target_sum_norm = unsafe { eval_degree_18(p_norm, r, w_i) };
    }
    // Final
    unsafe {
        let final_eq = eval_final_eq(&proof.r_x, &proof.r_y, tau_0);
        let final_alpha = eval_final_alpha(&proof.r_y, alpha);
        let pos_bound = 5 * (1 << 13) * 1024;
        let bal_bound = 9 * (1 << 13) * 1024;
        let final_ind_pos = eval_step_indicator(pos_bound, &proof.r_x, &proof.r_y);
        
        let ind_to_9 = eval_step_indicator(bal_bound, &proof.r_x, &proof.r_y);
        let final_ind_bal = ext_sub(&ind_to_9, &final_ind_pos);
        let final_w = proof.final_w;
        let final_m = proof.final_m;
        let mut ind_pos_8 = [0u32; 4];
        for k in 0..4 {
            let val = final_ind_pos[k] as u64 * 8;
            let mut res = (val as u32 as u64) + (val >> 32) * 99;
            if res >= q as u64 { res -= q as u64; }
            ind_pos_8[k] = res as u32;
        }
        let final_w_orig = extadd(&final_w, &ind_pos_8);
        let final_eval_alg = ext_mul_val(&ext_mul_val(&final_w_orig, &final_alpha), &final_m);

        let mut w_sq = ext_mul_val(&final_w, &final_w);
        let mut w_plus_8 = final_w;
        w_plus_8[0] = if w_plus_8[0] + 8 >= q { w_plus_8[0] + 8 - q } else { w_plus_8[0] + 8 };
        let mut c_bal = ext_mul_val(&final_w, &w_plus_8);
        
        for &k_sq in &[1, 4, 9, 16, 25, 36, 49] {
            let mut tmp = w_sq;
            tmp[0] = if tmp[0] >= k_sq { tmp[0] - k_sq } else { tmp[0] + q - k_sq };
            c_bal = ext_mul_val(&c_bal, &tmp);
        }
        
        let sum_ind = extadd(&final_ind_pos, &final_ind_bal);
        let final_eval_norm = ext_mul_val(&final_eq, &ext_mul_val(&c_bal, &sum_ind));
        if final_eval_alg == target_sum_alg && final_eval_norm == target_sum_norm {
            true
        } else {
            false
        }
    }
}