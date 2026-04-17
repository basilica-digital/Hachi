use crate::math::fields::*;
use crate::math::eq::build_eq_tbl_split;

pub unsafe fn eval_final_eq(r_x: &[[u32; 4]], r_y: &[[u32; 4]], tau_0: &[[u32; 4]]) -> [u32; 4] {
    let (mut eq_low, mut eq_high) = build_eq_tbl_split(tau_0);
    let ext_one = [1u32, 0, 0, 0];
    let mut cumulative_high = [1u32, 0, 0, 0];
    let mut h_len = 1 << 18;
    for round in 0..18 {
        let r = &r_x[round];
        let h_half = h_len / 2;
        
        let mut table_sum = [0u32; 4];
        for i in 0..h_half {
            eq_high[i] = fold_val(&eq_high[i], &eq_high[h_half + i], r);
            table_sum = extadd(&table_sum, &eq_high[i]);
        }
        h_len = h_half;
        
        let mut found_idx = -1;
        for test_idx in 0..28 {
            let t_val = &tau_0[test_idx];
            let r_t = ext_mul_val(r, t_val);
            let mut factor = extadd(&ext_one, &extadd(&r_t, &r_t));
            factor = ext_sub(&factor, r);
            factor = ext_sub(&factor, t_val);
            
            let test_cumulative = ext_mul_val(&cumulative_high, &factor);
            if test_cumulative == table_sum {
                found_idx = test_idx as i32;
                cumulative_high = test_cumulative;
                break;
            }
        }
    }

    let mut cumulative_low = [1u32, 0, 0, 0];
    let mut l_len = 1 << 10;
    for round in 0..10 {
        let r = &r_y[round];
        let l_half = l_len / 2;
        
        let mut table_sum = [0u32; 4];
        for i in 0..l_half {
            eq_low[i] = fold_val(&eq_low[i], &eq_low[l_half + i], r);
            table_sum = extadd(&table_sum, &eq_low[i]);
        }
        l_len = l_half;
        
        let mut found_idx = -1;
        for test_idx in 0..28 {
            let t_val = &tau_0[test_idx];
            let r_t = ext_mul_val(r, t_val);
            let mut factor = extadd(&ext_one, &extadd(&r_t, &r_t));
            factor = ext_sub(&factor, r);
            factor = ext_sub(&factor, t_val);
            
            let test_cumulative = ext_mul_val(&cumulative_low, &factor);
            if test_cumulative == table_sum {
                found_idx = test_idx as i32;
                cumulative_low = test_cumulative;
                break;
            }
        }
    }
    let final_res = ext_mul_val(&eq_high[0], &eq_low[0]);
    
    final_res
}

pub fn eval_final_alpha(r_y: &[[u32; 4]], alpha: &[u32; 4]) -> [u32; 4] {
    let mut res = [1u32, 0, 0, 0];
    let ext_one = [1u32, 0, 0, 0];

    let mut pows = [[0u32; 4]; 10];
    let mut curr = *alpha;
    for k in 0..10 {
        pows[k] = curr;
        curr = ext_mul_val(&curr, &curr); // alpha^{2^k}
    }

    for i in 0..10 {
        let r_i = &r_y[i];
        let a_pow = &pows[9 - i];

        // factor = (1 - r_i) + r_i * alpha^{2^k}
        let term1 = ext_sub(&ext_one, r_i);
        let term2 = ext_mul_val(r_i, a_pow);
        let factor = extadd(&term1, &term2);

        res = ext_mul_val(&res, &factor);
    }
    res
}

pub fn eval_step_indicator(k_bound: usize, r_x: &[[u32; 4]], r_y: &[[u32; 4]]) -> [u32; 4] {
    let mut sum = [0u32, 0, 0, 0];
    let mut prefix_eq = [1u32, 0, 0, 0];
    let ext_one = [1u32, 0, 0, 0];

    for j in 0..28 {
        let r_j = if j < 18 { &r_x[j] } else { &r_y[j - 18] };
        let bit = (k_bound >> (27 - j)) & 1;

        let one_minus_r = ext_sub(&ext_one, r_j);

        if bit == 1 {
            let branch_val = ext_mul_val(&prefix_eq, &one_minus_r);
            sum = extadd(&sum, &branch_val);
            prefix_eq = ext_mul_val(&prefix_eq, r_j);
        } else {
            prefix_eq = ext_mul_val(&prefix_eq, &one_minus_r);
        }
    }
    sum
}