use std::arch::x86_64::*;
use crate::math::fields::*;
use crate::math::field_simd::*;
use crate::mlp::eval::*;
use crate::Q;
use crate::hachi::prove::SumcheckProof;
use crate::hachi::setup::SetupParams;
use crate::hachi::verify::Challenge;
use crate::mlp::mle::{MBundle, mle_m};
use crate::math::eq::build_eq_tbl;
use std::time::Instant;

pub unsafe fn sumcheck_verify(
    proof: &SumcheckProof,
    params: &SetupParams,
    mut target_sum_alg: [u32; 4],
    mut target_sum_norm: [u32; 4],
    w_i: &[u32; 7],
    challenge: &Challenge
) -> bool {
    let q = 4294967197u32;
    let a0 = &challenge.a0;
    let a1 = &challenge.a1;
    let a2 = &challenge.a2;
    let a3 = &challenge.a3;
    let b0 = &challenge.b0;
    let b1 = &challenge.b1;
    let b2 = &challenge.b2;
    let b3 = &challenge.b3;
    let c = &challenge.c;
    let alpha = &challenge.alpha;
    let tau_0 = &challenge.tau_0;
    let tau_1 = &challenge.tau_1;

    let start = Instant::now();
    
    // X
    for round in 0..18 {
        let p_alg = &proof.p_alg_x[round];
        let p_norm = &proof.p_norm_x[round];
        let r = &proof.r_x[round]; 

        let current_sum_alg = unsafe { extadd(&p_alg[0], &p_alg[1]) };
        if current_sum_alg != target_sum_alg { return false; }
        
        let current_sum_norm = unsafe { extadd(&p_norm[0], &p_norm[1]) };
        if current_sum_norm != target_sum_norm { return false; }
        
        target_sum_alg = unsafe { eval_quadratic(&p_alg[0], &p_alg[1], &p_alg[2], r) };
        target_sum_norm = unsafe { eval_degree_6(p_norm, r, w_i) };

        
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
        target_sum_norm = unsafe { eval_degree_6(p_norm, r, w_i) };
    }

    // Final
    unsafe {
        let final_eq = eval_final_eq(&proof.r_x, &proof.r_y, tau_0);
        let final_alpha = eval_final_alpha(&proof.r_y, alpha[0..4].try_into().unwrap());
        let pos_bound = 5 * (1 << 14) * 1024;
        let bal_bound = 9 * (1 << 14) * 1024;
        let k_bound   = bal_bound + 64 * 1024;

        let ind_to_5 = eval_step_indicator(pos_bound, &proof.r_x, &proof.r_y);
        let ind_to_9 = eval_step_indicator(bal_bound, &proof.r_x, &proof.r_y);
        let ind_to_k = eval_step_indicator(k_bound, &proof.r_x, &proof.r_y);

        let final_ind_bal = ext_sub(&ind_to_9, &ind_to_5);
        let ind_k = ext_sub(&ind_to_k, &ind_to_9);
        let final_ind_pos = extadd(&ind_to_5, &ind_k);

        let final_w = proof.final_w;
        let height_2_ext = params.height_2*4;
        let height_z_ext = params.height_4*8*4;

        // M(alpha)
        let mut D_alpha = vec![0u32; height_2_ext];
        let mut E0_alpha = vec![0u32; 64];
        let mut E1_alpha = vec![0u32; 64];
        let mut E2_alpha = vec![0u32; 64];
        let mut E3_alpha = vec![0u32; 64];
        let mut DJ_alpha = vec![0u32; height_z_ext];

        let mut a0GJ_alpha = vec![0u32; height_z_ext];
        let mut a1GJ_alpha = vec![0u32; height_z_ext];
        let mut a2GJ_alpha = vec![0u32; height_z_ext];
        let mut a3GJ_alpha = vec![0u32; height_z_ext];
        let mut b0G_alpha = vec![0u32; height_2_ext];
        let mut b1G_alpha = vec![0u32; height_2_ext];
        let mut b2G_alpha = vec![0u32; height_2_ext];
        let mut b3G_alpha = vec![0u32; height_2_ext];
        let mut cG_alpha = vec![0u32; height_2_ext];

        let mut alpha_vec = vec![0u32; 4*1025];
        alpha_vec[0] = 1;
        alpha_vec[1] = 0;
        alpha_vec[2] = 0;
        alpha_vec[3] = 0;
        for i in 1..1025{
            let (prev, current) = alpha_vec.split_at_mut(4*i);
            extmul(&mut current[0..4], &prev[4*(i-1)..4*i], alpha);
        }
        // alpha_vec[4*1024] = - (alpha^1024 + 1)
        alpha_vec[4096] = ((params.q64 - ((alpha_vec[4096] as u64 + 1u64) % params.q64)) % params.q64) as u32;
        alpha_vec[4097] = ((params.q64 - (alpha_vec[4097] as u64 % params.q64)) % params.q64) as u32;
        alpha_vec[4098] = ((params.q64 - (alpha_vec[4098] as u64 % params.q64)) % params.q64) as u32;
        alpha_vec[4099] = ((params.q64 - (alpha_vec[4099] as u64 % params.q64)) % params.q64) as u32;

        let mut alpha_soa = [
            vec![0u32; params.n],
            vec![0u32; params.n],
            vec![0u32; params.n],
            vec![0u32; params.n],
        ];

        for i in 0..params.n {
            alpha_soa[0][i] = alpha_vec[4 * i];
            alpha_soa[1][i] = alpha_vec[4 * i + 1];
            alpha_soa[2][i] = alpha_vec[4 * i + 2];
            alpha_soa[3][i] = alpha_vec[4 * i + 3];
        }
        
        let tmp = alpha_vec[0..4].to_vec(); 
        
        let v_99 = _mm512_set1_epi64(99);
        let v_mask32 = _mm512_set1_epi64(0xFFFFFFFF);

        let load_ext = |ptr: *const u32| -> __m512i {
            _mm512_cvtepu32_epi64(_mm256_loadu_si256(ptr as *const __m256i))
        };

        let collapse_and_reduce = |simd_val: __m512i| -> u32 {
                let mut arr = [0u64; 8];
                _mm512_storeu_si512(arr.as_mut_ptr() as *mut __m512i, simd_val);
                let mut total_sum = 0u64;
                for x in arr { total_sum += x; }
                final_reduce(total_sum)
            };

        for i in 0..params.height_2 {
            let row_offset = i * params.n;
            let mut sum_d = [_mm512_setzero_si512(); 4];

            for j in (0..params.n).step_by(8) {
                let d_v = load_ext(params.d.as_ptr().add(row_offset + j));
                for c in 0..4 {
                    let a_c = load_ext(alpha_soa[c].as_ptr().add(j));
                    sum_d[c] = basemul_vec_register(sum_d[c], a_c, d_v, v_99, v_mask32);
                }
            }
            for c in 0..4 {
                D_alpha[i*4+c] = collapse_and_reduce(sum_d[c]);
            }
        }

        for i in 0..params.height_4 {
            let d_i = [D_alpha[4*i], D_alpha[4*i+1], D_alpha[4*i+2], D_alpha[4*i+3]];
            for j in 0..8 {
                let idx = i * 8 + j;
                let gv = params.g_matrix_2[j] as u64;
                for c in 0..4 {
                    let val = ((d_i[c] as u64 * gv) % params.q64) as u32;
                    DJ_alpha[4*idx + c] = if val == 0 { 0 } else { params.q64 as u32 - val };
                }
            }
        }

        for i in 0..16 {
            let off0 = i * params.n;          let off1 = (16 + i) * params.n;
            let off2 = (32 + i) * params.n;   let off3 = (48 + i) * params.n;

            let mut sum_e0 = [_mm512_setzero_si512(); 4]; let mut sum_e1 = [_mm512_setzero_si512(); 4];
            let mut sum_e2 = [_mm512_setzero_si512(); 4]; let mut sum_e3 = [_mm512_setzero_si512(); 4];

            for j in (0..params.n).step_by(8) {
                let e0_v = load_ext(params.e.as_ptr().add(off0 + j)); let e1_v = load_ext(params.e.as_ptr().add(off1 + j));
                let e2_v = load_ext(params.e.as_ptr().add(off2 + j)); let e3_v = load_ext(params.e.as_ptr().add(off3 + j));

                for c in 0..4 {
                    let a_c = load_ext(alpha_soa[c].as_ptr().add(j));
                    sum_e0[c] = basemul_vec_register(sum_e0[c], a_c, e0_v, v_99, v_mask32);
                    sum_e1[c] = basemul_vec_register(sum_e1[c], a_c, e1_v, v_99, v_mask32);
                    sum_e2[c] = basemul_vec_register(sum_e2[c], a_c, e2_v, v_99, v_mask32);
                    sum_e3[c] = basemul_vec_register(sum_e3[c], a_c, e3_v, v_99, v_mask32);
                }
            }
            for c in 0..4 {
                E0_alpha[i*4+c] = collapse_and_reduce(sum_e0[c]); E1_alpha[i*4+c] = collapse_and_reduce(sum_e1[c]);
                E2_alpha[i*4+c] = collapse_and_reduce(sum_e2[c]); E3_alpha[i*4+c] = collapse_and_reduce(sum_e3[c]);
            }
        }

        for i in 0..params.n{
            for j in 0..16{
                let gv = params.g_matrix_2[j];
                let ngv = ((params.q64 - gv as u64) % params.q64) as u32;
                for k in 0..16{
                    if c[i*16+k] < 1024{
                        let idx = c[i*16+k] as usize;
                        ext_base_mla(&mut cG_alpha[4*(i*16+j)..4*(i*16+j+1)], &alpha_vec[4*idx..4*(idx+1)], gv);
                    }else{
                        let idx = (c[i*16+k] - 1024) as usize;
                        ext_base_mla(&mut cG_alpha[4*(i*16+j)..4*(i*16+j+1)], &alpha_vec[4*idx..4*(idx+1)], ngv);
                    }
                }
            }
        }

        for i in 0..params.n {
            for j in 0..16 {
                let idx = i * 16 + j;
                let gv = params.g_matrix_2[j] as u64;
                b0G_alpha[4 * idx] = ((b0[i] as u64 * gv) % params.q64) as u32;
                b1G_alpha[4 * idx] = ((b1[i] as u64 * gv) % params.q64) as u32;
                b2G_alpha[4 * idx] = ((b2[i] as u64 * gv) % params.q64) as u32;
                b3G_alpha[4 * idx] = ((b3[i] as u64 * gv) % params.q64) as u32;
            }
        }
        
        for i in 0..params.n {
            let a0_u = a0[i] as u64;
            let a1_u = a1[i] as u64;
            let a2_u = a2[i] as u64;
            let a3_u = a3[i] as u64;
            let idx0 = i * 64;
            
            for j in 0..8 {
                let gv1 = params.g_matrix_4[j] as u64;
                let idx1 = idx0 + j * 8;
                for k in 0..8 {
                    let gv2 = (gv1 * params.g_matrix_2[k] as u64) % params.q64;
                    let idx2 = idx1 + k;
                    a0GJ_alpha[4 * idx2] = ((a0_u * gv2) % params.q64) as u32;
                    a1GJ_alpha[4 * idx2] = ((a1_u * gv2) % params.q64) as u32;
                    a2GJ_alpha[4 * idx2] = ((a2_u * gv2) % params.q64) as u32;
                    a3GJ_alpha[4 * idx2] = ((a3_u * gv2) % params.q64) as u32;
                }
            }
        }

        let tbl_tau1 = unsafe { build_eq_tbl(&challenge.tau_1) };

        let final_m = unsafe {
            let build_eq_msb_first = |r_slice: &[[u32; 4]]| -> Vec<[u32; 4]> {
                let mut tbl = vec![[1u32, 0, 0, 0]];
                let ext_one = [1u32, 0, 0, 0];
                for r_i in r_slice {
                    let one_minus_r = ext_sub(&ext_one, r_i);
                    let mut next_tbl = Vec::with_capacity(tbl.len() * 2);
                    for val in &tbl {
                        next_tbl.push(ext_mul_val(val, &one_minus_r));
                        next_tbl.push(ext_mul_val(val, r_i));
                    }
                    tbl = next_tbl;
                }
                tbl
            };

            let eq_4 = build_eq_msb_first(&proof.r_x[0..4]);
            let eq_14 = build_eq_msb_first(&proof.r_x[4..18]);

            let extdbl = |x: [u32; 4]| -> [u32; 4] { extadd(&x, &x) };
            let extneg = |x: [u32; 4]| -> [u32; 4] { ext_sub(&[0u32, 0, 0, 0], &x) };
            let w_tau_eq = |i: usize, u_top: usize| ext_mul_val(&tbl_tau1[i], &eq_4[u_top]);
            let add_all = |terms: &[[u32; 4]]| -> [u32; 4] {
                let mut sum = [0u32; 4];
                for t in terms { sum = extadd(&sum, t); }
                sum
            };
            
            let w_d = add_all(&[w_tau_eq(0, 0), w_tau_eq(1, 1), w_tau_eq(2, 2), w_tau_eq(3, 3), w_tau_eq(5, 4)]);
            let w_cg = add_all(&[w_tau_eq(10, 0), w_tau_eq(11, 1), w_tau_eq(12, 2), w_tau_eq(13, 3), w_tau_eq(14, 4)]);
            
            let w_b0g = add_all(&[w_tau_eq(6, 0), w_tau_eq(7, 1), w_tau_eq(8, 2), w_tau_eq(9, 3)]);
            let w_b1g = add_all(&[extdbl(w_tau_eq(6, 3)), w_tau_eq(7, 0), w_tau_eq(8, 1), w_tau_eq(9, 2)]);
            let w_b2g = add_all(&[extdbl(w_tau_eq(6, 2)), extdbl(w_tau_eq(7, 3)), w_tau_eq(8, 0), w_tau_eq(9, 1)]);
            let w_b3g = add_all(&[extdbl(w_tau_eq(6, 1)), extdbl(w_tau_eq(7, 2)), extdbl(w_tau_eq(8, 3)), w_tau_eq(9, 0)]);

            let mut w_a0gj = [[0u32; 4]; 4]; let mut w_a1gj = [[0u32; 4]; 4];
            let mut w_a2gj = [[0u32; 4]; 4]; let mut w_a3gj = [[0u32; 4]; 4];
            let mut w_dj   = [[0u32; 4]; 4];
            for c in 0..4 {
                w_a0gj[c] = extneg(w_tau_eq(10, 5 + c));
                w_a1gj[c] = extneg(w_tau_eq(11, 5 + c));
                w_a2gj[c] = extneg(w_tau_eq(12, 5 + c));
                w_a3gj[c] = extneg(w_tau_eq(13, 5 + c));
                w_dj[c]   = w_tau_eq(14, 5 + c);
            }
            
            let mut m_accum = [0u32; 4];

            let dot_product_simd = |arr: &Vec<u32>, offset: usize| -> [u32; 4] {
                let base = offset * 16384;
                let mut acc_soa: EF8 = (
                    _mm512_setzero_si512(), 
                    _mm512_setzero_si512(), 
                    _mm512_setzero_si512(), 
                    _mm512_setzero_si512()
                );
                
                let mut i = 0;
                while i < 16384 {
                    let arr_ptr = arr.as_ptr().add((base + i) * 4) as *const [u32; 4];
                    let eq_ptr = eq_14.as_ptr().add(i);
                    let v_arr = load8_soa(arr_ptr);
                    let v_eq = load8_soa(eq_ptr);
                    let v_mul = extmul_8x_local(v_arr, v_eq);
                    acc_soa = extadd_8x(acc_soa, v_mul);
                    
                    i += 8;
                }
                sum_soa_to_aos(acc_soa)
            };

            let eval_d = dot_product_simd(&D_alpha, 0);
            let eval_cg = dot_product_simd(&cG_alpha, 0);
            let eval_b0g = dot_product_simd(&b0G_alpha, 0);
            let eval_b1g = dot_product_simd(&b1G_alpha, 0);
            let eval_b2g = dot_product_simd(&b2G_alpha, 0);
            let eval_b3g = dot_product_simd(&b3G_alpha, 0);
            
            m_accum = extadd(&m_accum, &ext_mul_val(&w_d, &eval_d));
            m_accum = extadd(&m_accum, &ext_mul_val(&w_cg, &eval_cg));
            m_accum = extadd(&m_accum, &ext_mul_val(&w_b0g, &eval_b0g));
            m_accum = extadd(&m_accum, &ext_mul_val(&w_b1g, &eval_b1g));
            m_accum = extadd(&m_accum, &ext_mul_val(&w_b2g, &eval_b2g));
            m_accum = extadd(&m_accum, &ext_mul_val(&w_b3g, &eval_b3g));

            for c in 0..4 {
                let eval_a0gj = dot_product_simd(&a0GJ_alpha, c);
                let eval_a1gj = dot_product_simd(&a1GJ_alpha, c);
                let eval_a2gj = dot_product_simd(&a2GJ_alpha, c);
                let eval_a3gj = dot_product_simd(&a3GJ_alpha, c);
                let eval_dj = dot_product_simd(&DJ_alpha, c);

                m_accum = extadd(&m_accum, &ext_mul_val(&w_a0gj[c], &eval_a0gj));
                m_accum = extadd(&m_accum, &ext_mul_val(&w_a1gj[c], &eval_a1gj));
                m_accum = extadd(&m_accum, &ext_mul_val(&w_a2gj[c], &eval_a2gj));
                m_accum = extadd(&m_accum, &ext_mul_val(&w_a3gj[c], &eval_a3gj));
                m_accum = extadd(&m_accum, &ext_mul_val(&w_dj[c], &eval_dj));
            }

            for i in 0..64 {
                let mut sum_64 = [0u32; 4];
                let gv = params.g_matrix_2[i % 16];
                let neg_g_ext = [if gv == 0 { 0 } else { 4294967197u32 - gv }, 0, 0, 0];

                if i < 16 {
                    let val: [u32; 4] = E0_alpha[4*i..4*(i+1)].try_into().unwrap();
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(4, 9), &val));
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(0, 9), &neg_g_ext));
                } else if i < 32 {
                    let val: [u32; 4] = E1_alpha[4*(i-16)..4*(i-16+1)].try_into().unwrap();
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(4, 9), &val));
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(1, 9), &neg_g_ext));
                } else if i < 48 {
                    let val: [u32; 4] = E2_alpha[4*(i-32)..4*(i-32+1)].try_into().unwrap();
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(4, 9), &val));
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(2, 9), &neg_g_ext));
                } else {
                    let val: [u32; 4] = E3_alpha[4*(i-48)..4*(i-48+1)].try_into().unwrap();
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(4, 9), &val));
                    sum_64 = extadd(&sum_64, &ext_mul_val(&w_tau_eq(3, 9), &neg_g_ext));
                }
                m_accum = extadd(&m_accum, &ext_mul_val(&sum_64, &eq_14[i]));
            }

            let alpha_term: [u32; 4] = alpha_vec[4096..4100].try_into().unwrap();
            let mut alpha_offset_sum = [0u32; 4];
            for i in 0..15 {
                let weight = ext_mul_val(&tbl_tau1[i], &eq_4[9]);
                let weighted_eq = ext_mul_val(&weight, &eq_14[64 + i]);
                alpha_offset_sum = extadd(&alpha_offset_sum, &weighted_eq);
            }
            m_accum = extadd(&m_accum, &ext_mul_val(&alpha_term, &alpha_offset_sum));

            m_accum
        };
        
        let ind_pos_2 = extadd(&final_ind_pos, &final_ind_pos);
        let final_w_orig = extadd(&final_w, &ind_pos_2);
        
        let final_eval_alg = ext_mul_val(&ext_mul_val(&final_w_orig, &final_alpha), &final_m);
        let w_sq = ext_mul_val(&final_w, &final_w);
        let s_poly = extadd(&w_sq, &final_w);
        let mut s_minus_2 = s_poly;
        s_minus_2[0] = if s_minus_2[0] >= 2 { s_minus_2[0] - 2 } else { s_minus_2[0] + q - 2 };
        
        let c_bal = ext_mul_val(&s_poly, &s_minus_2);
        
        let sum_ind = extadd(&final_ind_pos, &final_ind_bal);
        let final_eval_norm = ext_mul_val(&final_eq, &ext_mul_val(&c_bal, &sum_ind));
        
        if final_eval_alg == target_sum_alg && final_eval_norm == target_sum_norm {
            true
        } else {
            false
        }
    }
}