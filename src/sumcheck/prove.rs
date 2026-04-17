use std::arch::x86_64::*;
use crate::math::fields::*;
use crate::math::field_simd::*;
use crate::utils::random::generate_random_q_element;
use crate::hachi::prove::SumcheckProof;
use std::time::Instant;



pub unsafe fn sumcheck_prove(
    W_compact: &[u32],
    mut M_table: Vec<[u32; 4]>, 
    mut alpha_vec: Vec<[u32; 4]>,
    mut eq_table_norm_low: Vec<[u32; 4]>,
    mut eq_table_norm_high: Vec<[u32; 4]>,
    ind_pos_compact: &[u8],
    ind_bal_compact: &[u8],
    tau_0: &[[u32; 4]; 28]
) -> SumcheckProof {
    
    let mut debug_eq_high = [1u32, 0, 0, 0];
    let mut debug_eq_low = [1u32, 0, 0, 0];
    let ext_one = [1u32, 0, 0, 0];
    let zero_soa: EF8 = (_mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512());
    let q = 4294967197u32;
    let q_u64 = q as u64;

    let mut W_table: Vec<[u32; 4]> = vec![[0u32; 4]; 1 << 27];
    let mut ind_pos_table: Vec<[u32; 4]> = vec![[0u32; 4]; 1 << 27];
    let mut ind_bal_table: Vec<[u32; 4]> = vec![[0u32; 4]; 1 << 27];

    let mut lut_c_sum = vec![[0u32; 7]; 256];
    let q_u64 = q as u64; 
    for w0 in 0..4 {
        for w1 in 0..4 {
            for ip0 in 0..2 {
                for ip1 in 0..2 {
                    for ib0 in 0..2 {
                        for ib1 in 0..2 {
                            let packed0 = (w0 << 2) | (ip0 << 1) | ib0;
                            let packed1 = (w1 << 2) | (ip1 << 1) | ib1;
                            let idx = (packed0 << 4) | packed1;
                            
                            let map_w = |w: u32, ip: u32| -> u64 {
                                if ip == 1 {
                                    if w >= 2 { (w - 2) as u64 } else { w as u64 + q_u64 - 2 }
                                } else { w as u64 }
                            };
                            let w0_m = map_w(w0, ip0);
                            let w1_m = map_w(w1, ip1);
                            
                            let w_diff = (w1_m + q_u64 - w0_m) % q_u64;
                            let ip_diff = (ip1 as u64 + q_u64 - ip0 as u64) % q_u64;
                            let ib_diff = (ib1 as u64 + q_u64 - ib0 as u64) % q_u64;
                            
                            for x in 0..7 {
                                let x_u64 = x as u64;
                                let w_x = (w0_m + x_u64 * w_diff) % q_u64;
                                let ip_x = (ip0 as u64 + x_u64 * ip_diff) % q_u64;
                                let ib_x = (ib0 as u64 + x_u64 * ib_diff) % q_u64;
                                
                                let w_sq = (w_x * w_x) % q_u64;
                                let s_poly = (w_sq + w_x) % q_u64;
                                let s_minus_2 = if s_poly >= 2 { s_poly - 2 } else { s_poly + q_u64 - 2 };
                                let c_poly = (s_poly * s_minus_2) % q_u64;
                                let sum_ind = (ip_x + ib_x) % q_u64;
                                
                                lut_c_sum[idx as usize][x] = ((c_poly * sum_ind) % q_u64) as u32;
                            }
                        }
                    }
                }
            }
        }
    }

    let mut proof = SumcheckProof {
        p_alg_x: Vec::with_capacity(18), p_norm_x: Vec::with_capacity(18), r_x: Vec::with_capacity(18),
        p_alg_y: Vec::with_capacity(10), p_norm_y: Vec::with_capacity(10), r_y: Vec::with_capacity(10),
        final_w: [0; 4]
    };

    let mut x_len = 1 << 18;
    let mut y_len = 1 << 10;
    let zero_soa: EF8 = (_mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512(), _mm512_setzero_si512());

    let get_eq = |idx: usize, eq_h: &[[u32; 4]], eq_l: &[[u32; 4]]| -> [u32; 4] {
        ext_mul_val(&eq_h[idx >> 10], &eq_l[idx & 0x3FF]) 
    };

    let mut packed_state = vec![0u8; 1 << 28];
    for i in 0..packed_state.len() {
        let w = (W_compact[i] & 3) as u8;
        let ip = ind_pos_compact[i];
        let ib = ind_bal_compact[i];
        packed_state[i] = (w << 2) | (ip << 1) | ib;
    }

    let mut x_len = 1 << 18;

    for round in 0..18 {
        let half_x = x_len / 2;
        let mut p_alg_soa = [zero_soa; 3];
        let mut p_norm_soa = [zero_soa; 7];

        // let start = Instant::now();
        
        for i in 0..half_x {
            let m0_soa = broadcast_aos_to_soa(&M_table[i]);
            let m1_soa = broadcast_aos_to_soa(&M_table[half_x + i]);
            let m_diff_soa = extsub_8x(m1_soa, m0_soa);
            let mut m_evals = [zero_soa; 3];
            let mut m_X = m0_soa;
            
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
            if round == 0 {
                let zero_reg = _mm512_setzero_si512();
                let q_vec = _mm512_set1_epi64(q as i64);
                let one_vec = _mm512_set1_epi64(1);
                let two_vec = _mm512_set1_epi64(2);
                let three_vec = _mm512_set1_epi64(3);
                let load_w_mapped = |idx: usize| -> EF8 {
                    let mut arr = [0u32; 8];
                    for lane in 0..8 {
                        let w_raw = W_compact[idx + lane]; 
                        let ip = ind_pos_compact[idx + lane];
                        arr[lane] = if ip == 1 {
                            if w_raw >= 2 { w_raw - 2 } else { w_raw + q - 2 }
                        } else { w_raw }; 
                    }
                    (_mm512_cvtepu32_epi64(_mm256_loadu_si256(arr.as_ptr() as *const __m256i)), zero_reg, zero_reg, zero_reg)
                };
                
                let load_ind = |ind_c: &[u8], idx: usize| -> EF8 {
                    let mut arr = [0u32; 8];
                    for lane in 0..8 { arr[lane] = ind_c[idx + lane] as u32; }
                    (_mm512_cvtepu32_epi64(_mm256_loadu_si256(arr.as_ptr() as *const __m256i)), zero_reg, zero_reg, zero_reg)
                };

                while j < y_len {
                    let idx0 = i * y_len + j;
                    let idx1 = (half_x + i) * y_len + j;
                    
                    let alpha_soa = load8_soa(alpha_vec.as_ptr().add(j));
                    let w0_raw_base = _mm512_cvtepu32_epi64(_mm256_loadu_si256(W_compact.as_ptr().add(idx0) as *const __m256i));
                    let w1_raw_base = _mm512_cvtepu32_epi64(_mm256_loadu_si256(W_compact.as_ptr().add(idx1) as *const __m256i));
                    
                    let mut w_orig_X = (w0_raw_base, zero_reg, zero_reg, zero_reg);
                    let w1_orig_ef8 = (w1_raw_base, zero_reg, zero_reg, zero_reg);
                    let w_orig_diff = extsub_8x(w1_orig_ef8, w_orig_X);
                    let ip0_base = _mm512_cvtepu8_epi64(_mm_loadl_epi64(ind_pos_compact.as_ptr().add(idx0) as *const __m128i));
                    let ip1_base = _mm512_cvtepu8_epi64(_mm_loadl_epi64(ind_pos_compact.as_ptr().add(idx1) as *const __m128i));
                    let ib0_base = _mm512_cvtepu8_epi64(_mm_loadl_epi64(ind_bal_compact.as_ptr().add(idx0) as *const __m128i));
                    let ib1_base = _mm512_cvtepu8_epi64(_mm_loadl_epi64(ind_bal_compact.as_ptr().add(idx1) as *const __m128i));

                    let ip0 = (ip0_base, zero_reg, zero_reg, zero_reg);
                    let ip1 = (ip1_base, zero_reg, zero_reg, zero_reg);
                    let ib0 = (ib0_base, zero_reg, zero_reg, zero_reg);
                    let ib1 = (ib1_base, zero_reg, zero_reg, zero_reg);

                    let ip0_times_2 = extadd_8x(ip0, ip0);
                    let ip1_times_2 = extadd_8x(ip1, ip1);
                    let w0_mapped = extsub_8x(w_orig_X, ip0_times_2);
                    let w1_mapped = extsub_8x(w1_orig_ef8, ip1_times_2);
                    
                    let mut w_mapped_X = w0_mapped;
                    let w_mapped_diff = extsub_8x(w1_mapped, w0_mapped);

                    let eq_h_soa = broadcast_aos_to_soa(&eq_table_norm_high[idx0 >> 10]);
                    let eq_h1_soa = broadcast_aos_to_soa(&eq_table_norm_high[idx1 >> 10]);
                    let eq_l_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx0 & 0x3FF));
                    let eq_l1_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx1 & 0x3FF));
                    
                    let eq_0 = extmul_8x_local(eq_h_soa, eq_l_soa);
                    let eq_1 = extmul_8x_local(eq_h1_soa, eq_l1_soa);
                    let eq_diff = extsub_8x(eq_1, eq_0);
                    let eq_2 = extadd_8x(eq_1, eq_diff);
                    let mut eq_X = eq_0;

                    // H, dH, ddH
                    let ip_diff_soa = extsub_8x(ip1, ip0);
                    let ib_diff_soa = extsub_8x(ib1, ib0);
                    let S_0 = extadd_8x(ip0, ib0);
                    let S_diff = extadd_8x(ip_diff_soa, ib_diff_soa);
                    let S_1 = extadd_8x(S_0, S_diff);
                    let S_2 = extadd_8x(S_1, S_diff);

                    let mut H = ext_mul_base_8x(eq_0, S_0.0);
                    let mut dH = extsub_8x(ext_mul_base_8x(eq_1, S_1.0), H);
                    let ddH = extsub_8x(extsub_8x(ext_mul_base_8x(eq_2, S_2.0), ext_mul_base_8x(eq_1, S_1.0)), dH);
                    
                    let two_ef8 = (two_vec, zero_reg, zero_reg, zero_reg);

                    for x in 0..7 {
                        if x < 3 {
                            let term = ext_mul_base_8x(alpha_soa, w_orig_X.0);
                            inner_w_alpha[x] = extadd_8x(inner_w_alpha[x], term);
                            w_orig_X = extadd_8x(w_orig_X, w_orig_diff);
                        }

                        let w_sq = ext_mul_base_8x(w_mapped_X, w_mapped_X.0);
                        let s_poly = extadd_8x(w_sq, w_mapped_X);
                        let s_minus_2 = extsub_8x(s_poly, two_ef8);
                        let c_poly = ext_mul_base_8x(s_poly, s_minus_2.0);
                        
                        p_norm_soa[x] = extadd_8x(p_norm_soa[x], ext_mul_base_8x(H, c_poly.0));

                        w_mapped_X = extadd_8x(w_mapped_X, w_mapped_diff);
                        eq_X = extadd_8x(eq_X, eq_diff);
                        H = extadd_8x(H, dH); 
                        dH = extadd_8x(dH, ddH);
                    }
                    j += 8;
                }
            } else { 
                while j < y_len { 
                    let idx0 = i * y_len + j;
                    let idx1 = (half_x + i) * y_len + j;
                    unsafe {
                        _mm_prefetch(W_table.as_ptr().add(idx0 + 16) as *const i8, _MM_HINT_T0);
                        _mm_prefetch(W_table.as_ptr().add(idx1 + 16) as *const i8, _MM_HINT_T0);
                    }
                    
                    let w0_soa = load8_soa(W_table.as_ptr().add(idx0));
                    let w1_soa = load8_soa(W_table.as_ptr().add(idx1));
                    let w_diff_soa = extsub_8x(w1_soa, w0_soa);

                    let ip0 = load8_soa(ind_pos_table.as_ptr().add(idx0));
                    let ip1 = load8_soa(ind_pos_table.as_ptr().add(idx1));
                    let ib0 = load8_soa(ind_bal_table.as_ptr().add(idx0));
                    let ib1 = load8_soa(ind_bal_table.as_ptr().add(idx1));
                    
                    let is_pos_zero = is_all_zeros_8x(ip0) && is_all_zeros_8x(ip1);
                    let is_bal_zero = is_all_zeros_8x(ib0) && is_all_zeros_8x(ib1);
                    let alpha_soa = load8_soa(alpha_vec.as_ptr().add(j));

                    if is_pos_zero && is_bal_zero {
                        let mut w_X_A = w0_soa;
                        for x in 0..3 {
                            let term = extmul_8x_local(w_X_A, alpha_soa);
                            inner_w_alpha[x] = extadd_8x(inner_w_alpha[x], term);
                            w_X_A = extadd_8x(w_X_A, w_diff_soa);
                        }
                    } else {
                        let ip_diff_soa = if is_pos_zero { zero_soa } else { extsub_8x(ip1, ip0) };
                        let ib_diff_soa = if is_bal_zero { zero_soa } else { extsub_8x(ib1, ib0) };
                        
                        let ip0_2 = extadd_8x(ip0, ip0); 
                        let mut w_orig_X = extadd_8x(w0_soa, ip0_2);
                        let w_orig_diff = extadd_8x(w_diff_soa, extadd_8x(ip_diff_soa, ip_diff_soa));

                        let eq_h_soa = broadcast_aos_to_soa(&eq_table_norm_high[idx0 >> 10]);
                        let eq_h1_soa = broadcast_aos_to_soa(&eq_table_norm_high[idx1 >> 10]);
                        let eq_l_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx0 & 0x3FF));
                        let eq_l1_soa = load8_soa(eq_table_norm_low.as_ptr().add(idx1 & 0x3FF));
                        
                        let eq_0 = extmul_8x_local(eq_h_soa, eq_l_soa);
                        let eq_1 = extmul_8x_local(eq_h1_soa, eq_l1_soa);
                        let eq_2 = extadd_8x(eq_1, extsub_8x(eq_1, eq_0));

                        let S_0 = extadd_8x(ip0, ib0);
                        let S_diff = extadd_8x(ip_diff_soa, ib_diff_soa);
                        let S_1 = extadd_8x(S_0, S_diff);
                        let S_2 = extadd_8x(S_1, S_diff);

                        let mut H = extmul_8x_local(eq_0, S_0);
                        let mut dH = extsub_8x(extmul_8x_local(eq_1, S_1), H);
                        let ddH = extsub_8x(extsub_8x(extmul_8x_local(eq_2, S_2), extmul_8x_local(eq_1, S_1)), dH);
                        
                        let mut w_mapped_X = w0_soa;

                        inner_w_alpha[0] = extadd_8x(inner_w_alpha[0], extmul_8x_local(w_orig_X, alpha_soa));
                        w_orig_X = extadd_8x(w_orig_X, w_orig_diff);

                        let w_sq_0 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_0 = extadd_8x(w_sq_0, w_mapped_X);
                        let c_poly_0 = extmul_8x_local(s_poly_0, sub_const_8x(s_poly_0, 2));
                        p_norm_soa[0] = extadd_8x(p_norm_soa[0], extmul_8x_local(H, c_poly_0));

                        w_mapped_X = extadd_8x(w_mapped_X, w_diff_soa);
                        H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);

                        inner_w_alpha[1] = extadd_8x(inner_w_alpha[1], extmul_8x_local(w_orig_X, alpha_soa));
                        w_orig_X = extadd_8x(w_orig_X, w_orig_diff);


                        let w_sq_1 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_1 = extadd_8x(w_sq_1, w_mapped_X);
                        let c_poly_1 = extmul_8x_local(s_poly_1, sub_const_8x(s_poly_1, 2));
                        p_norm_soa[1] = extadd_8x(p_norm_soa[1], extmul_8x_local(H, c_poly_1));

                        w_mapped_X = extadd_8x(w_mapped_X, w_diff_soa);
                        H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);

                        inner_w_alpha[2] = extadd_8x(inner_w_alpha[2], extmul_8x_local(w_orig_X, alpha_soa));

                        let w_sq_2 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_2 = extadd_8x(w_sq_2, w_mapped_X);
                        let c_poly_2 = extmul_8x_local(s_poly_2, sub_const_8x(s_poly_2, 2));
                        p_norm_soa[2] = extadd_8x(p_norm_soa[2], extmul_8x_local(H, c_poly_2));

                        w_mapped_X = extadd_8x(w_mapped_X, w_diff_soa);
                        H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);

                        let w_sq_3 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_3 = extadd_8x(w_sq_3, w_mapped_X);
                        let c_poly_3 = extmul_8x_local(s_poly_3, sub_const_8x(s_poly_3, 2));
                        p_norm_soa[3] = extadd_8x(p_norm_soa[3], extmul_8x_local(H, c_poly_3));

                        w_mapped_X = extadd_8x(w_mapped_X, w_diff_soa);
                        H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);

                        let w_sq_4 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_4 = extadd_8x(w_sq_4, w_mapped_X);
                        let c_poly_4 = extmul_8x_local(s_poly_4, sub_const_8x(s_poly_4, 2));
                        p_norm_soa[4] = extadd_8x(p_norm_soa[4], extmul_8x_local(H, c_poly_4));

                        w_mapped_X = extadd_8x(w_mapped_X, w_diff_soa);
                        H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);

                        let w_sq_5 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_5 = extadd_8x(w_sq_5, w_mapped_X);
                        let c_poly_5 = extmul_8x_local(s_poly_5, sub_const_8x(s_poly_5, 2));
                        p_norm_soa[5] = extadd_8x(p_norm_soa[5], extmul_8x_local(H, c_poly_5));

                        w_mapped_X = extadd_8x(w_mapped_X, w_diff_soa);
                        H = extadd_8x(H, dH); dH = extadd_8x(dH, ddH);

                        let w_sq_6 = extmul_8x_local(w_mapped_X, w_mapped_X);
                        let s_poly_6 = extadd_8x(w_sq_6, w_mapped_X);
                        let c_poly_6 = extmul_8x_local(s_poly_6, sub_const_8x(s_poly_6, 2));
                        p_norm_soa[6] = extadd_8x(p_norm_soa[6], extmul_8x_local(H, c_poly_6));
                    }
                    j += 8;
                }
            }
            for x in 0..3 {
                p_alg_soa[x] = extadd_8x(p_alg_soa[x], extmul_8x_local(inner_w_alpha[x], m_evals[x]));
            }
        }

        // let duration = start.elapsed();
        // println!("- sumcheck evaluation: {:?}", duration);
        // let start = Instant::now();

        let mut p_alg = [[0u32; 4]; 3]; let mut p_norm = [[0u32; 4]; 7];
        for x in 0..3 { p_alg[x] = sum_soa_to_aos(p_alg_soa[x]); }
        for x in 0..7 { p_norm[x] = sum_soa_to_aos(p_norm_soa[x]); }
        
        proof.p_alg_x.push(p_alg);
        proof.p_norm_x.push(p_norm);
        
        let r_raw = generate_random_q_element(1, 16);
        let mut r = [0u32; 4]; r.copy_from_slice(&r_raw[0..4]);
        proof.r_x.push(r);
        let r_soa = broadcast_aos_to_soa(&r);
        let zero_reg = _mm512_setzero_si512();

        for i in 0..half_x {
            M_table[i] = fold_val(&M_table[i], &M_table[half_x + i], &r);
            eq_table_norm_high[i] = fold_val(&eq_table_norm_high[i], &eq_table_norm_high[half_x + i], &r);
            let mut j = 0;
            if round == 0 {
                let zero_reg = _mm512_setzero_si512();
                let q_vec = _mm512_set1_epi64(q as i64);
                let two_vec = _mm512_set1_epi64(2);
                let one_vec = _mm512_set1_epi64(1);

                while j < y_len {
                    let idx0 = i * y_len + j; 
                    let idx1 = (half_x + i) * y_len + j;
                    
                    let load_w_mapped = |idx: usize| -> EF8 {
                        let w_raw = _mm512_cvtepu32_epi64(_mm256_loadu_si256(W_compact.as_ptr().add(idx) as *const __m256i));
                        let ip_128 = _mm_loadl_epi64(ind_pos_compact.as_ptr().add(idx) as *const __m128i);
                        let ip_512 = _mm512_cvtepu8_epi64(ip_128);
                        
                        let mask_ip_1 = _mm512_cmpeq_epi64_mask(ip_512, one_vec);
                        let mask_ge_2 = _mm512_cmpge_epu64_mask(w_raw, two_vec);
                        
                        let w_minus_2 = _mm512_sub_epi64(w_raw, two_vec);
                        let w_plus_q_minus_2 = _mm512_add_epi64(w_minus_2, q_vec);
                        
                        let w_ip_1 = _mm512_mask_blend_epi64(mask_ge_2, w_plus_q_minus_2, w_minus_2);
                        let w_mapped = _mm512_mask_blend_epi64(mask_ip_1, w_raw, w_ip_1);
                        
                        (w_mapped, zero_reg, zero_reg, zero_reg)
                    };
                    
                    let load_ind = |ind_c: &[u8], idx: usize| -> EF8 {
                        let ind_128 = _mm_loadl_epi64(ind_c.as_ptr().add(idx) as *const __m128i);
                        (_mm512_cvtepu8_epi64(ind_128), zero_reg, zero_reg, zero_reg)
                    };

                    store8_soa(W_table.as_mut_ptr().add(idx0), fold_8x(load_w_mapped(idx0), load_w_mapped(idx1), r_soa));
                    store8_soa(ind_pos_table.as_mut_ptr().add(idx0), fold_8x(load_ind(ind_pos_compact, idx0), load_ind(ind_pos_compact, idx1), r_soa));
                    store8_soa(ind_bal_table.as_mut_ptr().add(idx0), fold_8x(load_ind(ind_bal_compact, idx0), load_ind(ind_bal_compact, idx1), r_soa));
                    j += 8;
                }
            } else {
                while j < y_len {
                    let idx0 = i * y_len + j; let idx1 = (half_x + i) * y_len + j;
                    store8_soa(W_table.as_mut_ptr().add(idx0), fold_8x(load8_soa(W_table.as_ptr().add(idx0)), load8_soa(W_table.as_ptr().add(idx1)), r_soa));
                    store8_soa(ind_pos_table.as_mut_ptr().add(idx0), fold_8x(load8_soa(ind_pos_table.as_ptr().add(idx0)), load8_soa(ind_pos_table.as_ptr().add(idx1)), r_soa));
                    store8_soa(ind_bal_table.as_mut_ptr().add(idx0), fold_8x(load8_soa(ind_bal_table.as_ptr().add(idx0)), load8_soa(ind_bal_table.as_ptr().add(idx1)), r_soa));
                    j += 8;
                }
            }
        }

        // let duration = start.elapsed();
        // println!("- sumcheck fold: {:?}", duration);
        // let start = Instant::now();

        x_len = half_x;
    }

    let m_scalar = M_table[0]; 
    for round in 0..10 {
        let half_y = y_len / 2;
        let mut p_alg = [[0u32; 4]; 3];
        let mut p_norm = [[0u32; 4]; 7];
        
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
            let mut alpha_m_X = ext_mul_val(&m_scalar, &a0); let mut alpha_m_diff = ext_mul_val(&m_scalar, &a_diff);
            let mut w_X_A = w0; let mut ind_pos_X_A = ind_pos_0;
            
            for x in 0..3 {
                let w_orig_X = if is_pos_zero {
                    w_X_A
                } else {
                    let ind_pos_2 = extadd(&ind_pos_X_A, &ind_pos_X_A);
                    extadd(&w_X_A, &ind_pos_2)
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
                
                for x in 0..7 {
                    let w_sq = ext_mul_val(&w_X, &w_X);
                    let s_poly = extadd(&w_sq, &w_X);
                    
                    let mut s_minus_2 = s_poly;
                    s_minus_2[0] = if s_minus_2[0] >= 2 { s_minus_2[0] - 2 } else { s_minus_2[0] + q - 2 };
                    
                    let c_poly = ext_mul_val(&s_poly, &s_minus_2);
                    let sum_ind = extadd(&ind_pos_X, &ind_bal_X);
                    
                    p_norm[x] = extadd(&p_norm[x], &ext_mul_val(&eq_X, &ext_mul_val(&c_poly, &sum_ind)));
                    
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
    proof
}