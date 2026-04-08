use std::arch::x86_64::*;
use crate::math::fields::*;
use crate::math::field_simd::*;
use crate::utils::random::generate_random_q_element;
use crate::hachi::prove::SumcheckProof;

pub unsafe fn sumcheck_prove(
    mut W_table: Vec<[u32; 4]>,
    mut M_table: Vec<[u32; 4]>, 
    mut alpha_vec: Vec<[u32; 4]>,
    mut eq_table_norm_low: Vec<[u32; 4]>,
    mut eq_table_norm_high: Vec<[u32; 4]>,
    mut ind_pos_table: Vec<[u32; 4]>,
    mut ind_bal_table: Vec<[u32; 4]>,
    tau_0: &[[u32; 4]; 27]
) -> SumcheckProof {
    
    let mut debug_eq_high = [1u32, 0, 0, 0];
    let mut debug_eq_low = [1u32, 0, 0, 0];
    let ext_one = [1u32, 0, 0, 0];

    let mut proof = SumcheckProof {
        p_alg_x: Vec::with_capacity(17), p_norm_x: Vec::with_capacity(17), r_x: Vec::with_capacity(17),
        p_alg_y: Vec::with_capacity(10), p_norm_y: Vec::with_capacity(10), r_y: Vec::with_capacity(10),
        final_w: [0; 4], final_alpha: [0; 4], final_m: [0; 4], final_eq: [0; 4], final_ind_pos: [0; 4], final_ind_bal: [0; 4],
    };

    let q = 4294967197u32;
    for i in 0..W_table.len() {
        if ind_pos_table[i] == [1, 0, 0, 0] {
            W_table[i][0] = if W_table[i][0] >= 8 { W_table[i][0] - 8 } else { W_table[i][0] + q - 8 };
        }
    }

    let mut x_len = 1 << 17;
    let mut y_len = 1 << 10;
    let zero_soa: EF8 = (_mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512());

    let get_eq = |idx: usize, eq_h: &[[u32; 4]], eq_l: &[[u32; 4]]| -> [u32; 4] {
        ext_mul_val(&eq_h[idx >> 10], &eq_l[idx & 0x3FF]) 
    };

    for round in 0..17 {
        let half_x = x_len / 2;
        let mut p_alg_soa = [zero_soa; 3];
        let mut p_norm_soa = [zero_soa; 19];
        
        for i in 0..half_x {
            let m0_soa = broadcast_aos_to_soa(&M_table[i]);
            let m1_soa = broadcast_aos_to_soa(&M_table[half_x + i]);
            let m_diff_soa = extsub_8x(m1_soa, m0_soa);
            let mut m_evals = [zero_soa; 3];
            let mut m_X = m0_soa;
            // M linear interpolation
            for x in 0..3 {
                m_evals[x] = m_X;
                m_X = extadd_8x(m_X, m_diff_soa);
            }

            let h_idx0 = (i * y_len) >> 10;
            let h_idx1 = ((half_x + i) * y_len) >> 10;
            let eq_h_soa = broadcast_aos_to_soa(&eq_table_norm_high[h_idx0]);
            let eq_h1_soa = broadcast_aos_to_soa(&eq_table_norm_high[h_idx1]);
            let mut inner_w_alpha = [zero_soa; 3];

            let mut j = 0;
            while j < y_len { 
                let idx0 = i * y_len + j;
                let idx1 = (half_x + i) * y_len + j;
                
                let w0_soa = load8_soa(W_table.as_ptr().add(idx0));
                let w1_soa = load8_soa(W_table.as_ptr().add(idx1));
                let w_diff_soa = extsub_8x(w1_soa, w0_soa);

                let ip0 = load8_soa(ind_pos_table.as_ptr().add(idx0));
                let ip1 = load8_soa(ind_pos_table.as_ptr().add(idx1));
                let ib0 = load8_soa(ind_bal_table.as_ptr().add(idx0));
                let ib1 = load8_soa(ind_bal_table.as_ptr().add(idx1));
                
                let is_pos_zero = is_all_zeros_8x(ip0) && is_all_zeros_8x(ip1);
                let is_bal_zero = is_all_zeros_8x(ib0) && is_all_zeros_8x(ib1);

                // Constraint
                let alpha_soa = load8_soa(alpha_vec.as_ptr().add(j));
                let mut w_X_A = w0_soa;
                let mut ip_X_A = ip0;
                let ip_diff_soa = if is_pos_zero { zero_soa } else { extsub_8x(ip1, ip0) };
                
                if is_pos_zero {
                    for x in 0..3 {
                        let term = if round == 0 { ext_mul_base_8x(alpha_soa, w_X_A.0) } else { extmul_8x_local(w_X_A, alpha_soa) };
                        inner_w_alpha[x] = extadd_8x(inner_w_alpha[x], term);
                        w_X_A = extadd_8x(w_X_A, w_diff_soa);
                    }
                } else {
                    for x in 0..3 {
                        let ip_8 = ext_scalar_mul_8_8x(ip_X_A);
                        let w_orig_X = extadd_8x(w_X_A, ip_8);
                        let term = if round == 0 { ext_mul_base_8x(alpha_soa, w_orig_X.0) } else { extmul_8x_local(w_orig_X, alpha_soa) };
                        inner_w_alpha[x] = extadd_8x(inner_w_alpha[x], term);
                        w_X_A = extadd_8x(w_X_A, w_diff_soa);
                        ip_X_A = extadd_8x(ip_X_A, ip_diff_soa);
                    }
                }

                // Norm
                if round == 0 {
                    let w0_base = w0_soa.0; let w_diff_base = w_diff_soa.0;
                    let ip0_base = ip0.0; let ip_diff_base = ip_diff_soa.0;
                    let ib0_base = ib0.0; let ib_diff_base = if is_bal_zero { _mm512_setzero_si512() } else { base_sub_8x(ib1.0, ib0.0) };

                    if !is_pos_zero || !is_bal_zero {
                        let eq_l_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx0 & 0x3FF));
                        let eq_0 = extmul_8x_local(eq_h_soa, eq_l_soa);

                        let eq_l1_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx1 & 0x3FF));
                        let eq_1 = extmul_8x_local(eq_h1_soa, eq_l1_soa);
                        let eq_diff_soa = extsub_8x(eq_1, eq_0);
                        let eq_2 = extadd_8x(eq_1, eq_diff_soa);

                        let S_0 = base_add_8x(ip0_base, ib0_base);
                        let S_diff = base_add_8x(ip_diff_base, ib_diff_base);
                        let S_1 = base_add_8x(S_0, S_diff);
                        let S_2 = base_add_8x(S_1, S_diff);

                        let H_0 = ext_mul_base_8x(eq_0, S_0);
                        let H_1 = ext_mul_base_8x(eq_1, S_1);
                        let H_2 = ext_mul_base_8x(eq_2, S_2);

                        let mut H = H_0;
                        let mut dH = extsub_8x(H_1, H_0);
                        let dH_1 = extsub_8x(H_2, H_1);
                        let ddH = extsub_8x(dH_1, dH);

                        let w_1 = base_add_8x(w0_base, w_diff_base);
                        let w_2 = base_add_8x(w_1, w_diff_base);

                        let Z_0 = base_mul_8x(w0_base, w0_base);
                        let Z_1 = base_mul_8x(w_1, w_1);
                        let Z_2 = base_mul_8x(w_2, w_2);

                        let mut Z = Z_0;
                        let mut dZ = base_sub_8x(Z_1, Z_0);
                        let dZ_1 = base_sub_8x(Z_2, Z_1);
                        let ddZ = base_sub_8x(dZ_1, dZ);

                        let mut w_8 = base_add_8x(base_add_8x(w0_base, w0_base), base_add_8x(w0_base, w0_base));
                        w_8 = base_add_8x(w_8, w_8); 
                        let dw_8 = base_add_8x(base_add_8x(w_diff_base, w_diff_base), base_add_8x(w_diff_base, w_diff_base));
                        let dw_8 = base_add_8x(dw_8, dw_8);

                        for x in 0..19 {
                            let Z2 = base_mul_8x(Z, Z);
                            let mut Q1 = base_sub_const_8x(Z, 50); Q1 = base_mul_8x(Q1, Z); Q1 = base_add_8x(Q1, _mm512_set1_epi64(49));
                            let mut Q2 = base_sub_const_8x(Z, 40); Q2 = base_mul_8x(Q2, Z); Q2 = base_add_8x(Q2, _mm512_set1_epi64(144));
                            let mut Q3 = base_sub_const_8x(Z, 34); Q3 = base_mul_8x(Q3, Z); Q3 = base_add_8x(Q3, _mm512_set1_epi64(225));
                            let Q4 = base_sub_const_8x(Z, 16);
                            
                            let Q12 = base_mul_8x(Q1, Q2);
                            let Q34 = base_mul_8x(Q3, Q4);
                            let P_Z = base_mul_8x(Q12, Q34);
                            
                            let c_bal_base = base_mul_8x(base_add_8x(Z, w_8), P_Z);
                            let term = ext_mul_base_8x(H, c_bal_base);
                            p_norm_soa[x] = extadd_8x(p_norm_soa[x], term);
                            
                            Z = base_add_8x(Z, dZ); dZ = base_add_8x(dZ, ddZ);
                            H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);
                            w_8 = base_add_8x(w_8, dw_8);
                        }
                    }
                } else {
                    if !is_pos_zero || !is_bal_zero {
                        let ib_diff_soa = extsub_8x(ib1, ib0);
                        let eq_l_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx0 & 0x3FF));
                        let eq_0 = extmul_8x_local(eq_h_soa, eq_l_soa);

                        let eq_l1_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx1 & 0x3FF));
                        let eq_1 = extmul_8x_local(eq_h1_soa, eq_l1_soa);
                        let eq_diff_soa = extsub_8x(eq_1, eq_0);
                        let eq_2 = extadd_8x(eq_1, eq_diff_soa);

                        let S_0 = extadd_8x(ip0, ib0);
                        let S_diff = extadd_8x(ip_diff_soa, ib_diff_soa);
                        let S_1 = extadd_8x(S_0, S_diff);
                        let S_2 = extadd_8x(S_1, S_diff);

                        // h = eq * ind
                        let H_0 = extmul_8x_local(eq_0, S_0);
                        let H_1 = extmul_8x_local(eq_1, S_1);
                        let H_2 = extmul_8x_local(eq_2, S_2);

                        let mut H = H_0;
                        let mut dH = extsub_8x(H_1, H_0);
                        let dH_1 = extsub_8x(H_2, H_1);
                        let ddH = extsub_8x(dH_1, dH);

                        let w_1 = extadd_8x(w0_soa, w_diff_soa);
                        let w_2 = extadd_8x(w_1, w_diff_soa);

                        let Z_0 = extmul_8x_local(w0_soa, w0_soa);
                        let Z_1 = extmul_8x_local(w_1, w_1);
                        let Z_2 = extmul_8x_local(w_2, w_2);

                        let mut Z = Z_0;
                        let mut dZ = extsub_8x(Z_1, Z_0);
                        let dZ_1 = extsub_8x(Z_2, Z_1);
                        let ddZ = extsub_8x(dZ_1, dZ);

                        let mut w_8 = ext_scalar_mul_8_8x(w0_soa);
                        let dw_8 = ext_scalar_mul_8_8x(w_diff_soa);
                        
                        for x in 0..19 {
                            let Z2 = extmul_8x_local(Z, Z); 
                            let Q1 = ext_fma_scalar_const_8x(Z2, Z, 4294967147, 49);  
                            let Q2 = ext_fma_scalar_const_8x(Z2, Z, 4294967157, 144); 
                            let Q3 = ext_fma_scalar_const_8x(Z2, Z, 4294967163, 225); 
                            let Q4 = sub_const_8x(Z, 16);
                            
                            let Q12 = extmul_8x_local(Q1, Q2);
                            let Q34 = extmul_8x_local(Q3, Q4);
                            let P_Z = extmul_8x_local(Q12, Q34);
                            
                            let c_bal = extmul_8x_local(extadd_8x(Z, w_8), P_Z);
                            let term = extmul_8x_local(H, c_bal);
                            
                            p_norm_soa[x] = extadd_8x(p_norm_soa[x], term);
                            
                            Z = extadd_8x(Z, dZ); dZ = extadd_8x(dZ, ddZ);
                            H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);
                            w_8 = extadd_8x(w_8, dw_8);
                        }
                    }
                }
                j += 8;
            }
            for x in 0..3 {
                p_alg_soa[x] = extadd_8x(p_alg_soa[x], extmul_8x_local(inner_w_alpha[x], m_evals[x]));
            }
        }

        let mut p_alg = [[0u32; 4]; 3]; let mut p_norm = [[0u32; 4]; 19];
        for x in 0..3 { p_alg[x] = sum_soa_to_aos(p_alg_soa[x]); }
        for x in 0..19 { p_norm[x] = sum_soa_to_aos(p_norm_soa[x]); }
        
        proof.p_alg_x.push(p_alg);
        proof.p_norm_x.push(p_norm);
        
        // r from transcript
        let r_raw = generate_random_q_element(1, 16);
        let mut r = [0u32; 4]; r.copy_from_slice(&r_raw[0..4]);
        proof.r_x.push(r);
        let r_soa = broadcast_aos_to_soa(&r);
        
        // folding
        for i in 0..half_x {
            M_table[i] = fold_val(&M_table[i], &M_table[half_x + i], &r);
            let mut j = 0;
            while j < y_len{
                let idx0 = i * y_len + j; let idx1 = (half_x + i) * y_len + j;
                store8_soa(W_table.as_mut_ptr().add(idx0), fold_8x(load8_soa(W_table.as_ptr().add(idx0)), load8_soa(W_table.as_ptr().add(idx1)), r_soa));
                store8_soa(ind_pos_table.as_mut_ptr().add(idx0), fold_8x(load8_soa(ind_pos_table.as_ptr().add(idx0)), load8_soa(ind_pos_table.as_ptr().add(idx1)), r_soa));
                store8_soa(ind_bal_table.as_mut_ptr().add(idx0), fold_8x(load8_soa(ind_bal_table.as_ptr().add(idx0)), load8_soa(ind_bal_table.as_ptr().add(idx1)), r_soa));
                j += 8;
            }
        }

        let bit_to_fold = 26 - round; 
        if bit_to_fold >= 10 {
            let h_half = 1 << (bit_to_fold - 10);
            for i in 0..h_half { eq_table_norm_high[i] = fold_val(&eq_table_norm_high[i], &eq_table_norm_high[h_half + i], &r); }
        } else {
            let l_half = 1 << bit_to_fold;
            for i in 0..l_half { eq_table_norm_low[i] = fold_val(&eq_table_norm_low[i], &eq_table_norm_low[l_half + i], &r); }
        }
        x_len = half_x;
    }

    let m_scalar = M_table[0]; 
    for round in 0..10 {
        let half_y = y_len / 2;
        let mut p_alg = [[0u32; 4]; 3];
        let mut p_norm = [[0u32; 4]; 19];
        
        for j in 0..half_y {
            let w0 = W_table[j]; let w1 = W_table[half_y + j]; 
            let w_diff = ext_sub(&w1, &w0);
            let ind_pos_0 = ind_pos_table[j]; let ind_pos_1 = ind_pos_table[half_y + j];
            let ind_bal_0 = ind_bal_table[j]; let ind_bal_1 = ind_bal_table[half_y + j];
            let zero_arr = [0u32; 4];
            let is_pos_zero = ind_pos_0 == zero_arr && ind_pos_1 == zero_arr;
            let is_bal_zero = ind_bal_0 == zero_arr && ind_bal_1 == zero_arr;
            let ind_pos_diff = if is_pos_zero { zero_arr } else { ext_sub(&ind_pos_1, &ind_pos_0) };
            
            let a0 = alpha_vec[j]; let a1 = alpha_vec[half_y + j]; let a_diff = ext_sub(&a1, &a0);
            let mut alpha_m_X = ext_mul_val(&m_scalar, &a0); let alpha_m_diff = ext_mul_val(&m_scalar, &a_diff);
            let mut w_X_A = w0; let mut ind_pos_X_A = ind_pos_0;
            
            for x in 0..3 {
                let w_orig_X = if is_pos_zero {
                    w_X_A
                } else {
                    let mut ind_pos_8 = [0u32; 4];
                    for k in 0..4 {
                        let val = ind_pos_X_A[k] as u64 * 8;
                        let mut res = (val as u32 as u64) + (val >> 32) * 99;
                        if res >= q as u64 { res -= q as u64; }
                        ind_pos_8[k] = res as u32;
                    }
                    extadd(&w_X_A, &ind_pos_8)
                };
                p_alg[x] = extadd(&p_alg[x], &ext_mul_val(&w_orig_X, &alpha_m_X));
                alpha_m_X = extadd(&alpha_m_X, &alpha_m_diff);
                w_X_A = extadd(&w_X_A, &w_diff);
                if !is_pos_zero { ind_pos_X_A = extadd(&ind_pos_X_A, &ind_pos_diff); }
            }
            
            if !is_pos_zero || !is_bal_zero {
                let eq0 = get_eq(j, &eq_table_norm_high, &eq_table_norm_low);
                let eq1 = get_eq(half_y + j, &eq_table_norm_high, &eq_table_norm_low);
                let eq_diff = ext_sub(&eq1, &eq0);
                let ind_bal_diff = if is_bal_zero { zero_arr } else { ext_sub(&ind_bal_1, &ind_bal_0) };
                let mut w_X = w0; let mut eq_X = eq0;
                let mut ind_pos_X = ind_pos_0; let mut ind_bal_X = ind_bal_0;
                
                for x in 0..19 {
                    let w_sq = ext_mul_val(&w_X, &w_X);
                    let mut w_8 = [0u32; 4];
                    for k in 0..4 { 
                        let val = w_X[k] as u64 * 8;
                        let mut res = (val as u32 as u64) + (val >> 32) * 99;
                        if res >= q as u64 { res -= q as u64; }
                        w_8[k] = res as u32;
                    }
                    let mut c_bal = extadd(&w_sq, &w_8);

                    for &k_sq in &[1, 4, 9, 16, 25, 36, 49] {
                        let mut w_sq_minus_k = w_sq;
                        w_sq_minus_k[0] = if w_sq_minus_k[0] >= k_sq { w_sq_minus_k[0] - k_sq } else { w_sq_minus_k[0] + q - k_sq };
                        c_bal = ext_mul_val(&c_bal, &w_sq_minus_k);
                    }
                    
                    let sum_ind = extadd(&ind_pos_X, &ind_bal_X);
                    p_norm[x] = extadd(&p_norm[x], &ext_mul_val(&eq_X, &ext_mul_val(&c_bal, &sum_ind)));
                    
                    w_X = extadd(&w_X, &w_diff); eq_X = extadd(&eq_X, &eq_diff);
                    if !is_pos_zero { ind_pos_X = extadd(&ind_pos_X, &ind_pos_diff); }
                    if !is_bal_zero { ind_bal_X = extadd(&ind_bal_X, &ind_bal_diff); }
                }
            }
        }

        proof.p_alg_y.push(p_alg);
        proof.p_norm_y.push(p_norm);
        let r_raw = generate_random_q_element(1, 16);
        let mut r = [0u32; 4]; r.copy_from_slice(&r_raw[0..4]);
        proof.r_y.push(r);
        
        for j in 0..half_y {
            W_table[j] = fold_val(&W_table[j], &W_table[half_y + j], &r);
            alpha_vec[j] = fold_val(&alpha_vec[j], &alpha_vec[half_y + j], &r);
            ind_pos_table[j] = fold_val(&ind_pos_table[j], &ind_pos_table[half_y + j], &r);
            ind_bal_table[j] = fold_val(&ind_bal_table[j], &ind_bal_table[half_y + j], &r);
        }

        let bit_to_fold = 9 - round;
        let l_half = 1 << bit_to_fold;
        for i in 0..l_half {
            eq_table_norm_low[i] = fold_val(&eq_table_norm_low[i], &eq_table_norm_low[l_half + i], &r);
        }
        y_len = half_y;
    }
    proof.final_w = W_table[0];
    proof.final_alpha = alpha_vec[0];
    proof.final_m = M_table[0];
    proof.final_eq = ext_mul_val(&eq_table_norm_high[0], &eq_table_norm_low[0]);
    proof.final_ind_pos = ind_pos_table[0];
    proof.final_ind_bal = ind_bal_table[0];
    proof
}