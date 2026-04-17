use std::arch::x86_64::*;
use rayon::prelude::*;
use crate::math::field_simd::extmul_8x;

// move
use crate::{ext_sub, extmul, u642u32};

unsafe fn build_eq_tbl_base(tau: &[[u32; 4]]) -> Vec<[u32; 4]> {
    let n = tau.len();
    let total_size = 1 << n;
    let mut eq_table = vec![[0u32; 4]; total_size];
    eq_table[0] = [1, 0, 0, 0];
    let ext_one = [1u32, 0, 0, 0];

    for j in 0..n {
        let current_len = 1 << j;
        let t_j = &tau[j];
        let one_minus_t_j = ext_sub(&ext_one, t_j);

        if current_len < 8 {
            for i in (0..current_len).rev() {
                let val = eq_table[i];
                let mut right = [0u32; 4];
                extmul(&mut right, &val, t_j);
                eq_table[i + current_len] = right;
                
                let mut left = [0u32; 4];
                extmul(&mut left, &val, &one_minus_t_j);
                eq_table[i] = left;
            }
        } else {
            let t_b0 = _mm512_set1_epi64(t_j[0] as i64);
            let t_b1 = _mm512_set1_epi64(t_j[1] as i64);
            let t_b2 = _mm512_set1_epi64(t_j[2] as i64);
            let t_b3 = _mm512_set1_epi64(t_j[3] as i64);

            let omt_b0 = _mm512_set1_epi64(one_minus_t_j[0] as i64);
            let omt_b1 = _mm512_set1_epi64(one_minus_t_j[1] as i64);
            let omt_b2 = _mm512_set1_epi64(one_minus_t_j[2] as i64);
            let omt_b3 = _mm512_set1_epi64(one_minus_t_j[3] as i64);

            if j < 14 {
                let mut i = current_len;
                while i >= 8 {
                    i -= 8;
                    let ptr = eq_table.as_ptr().add(i);
                    let v0 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[0] as i32, (*ptr.add(6))[0] as i32, (*ptr.add(5))[0] as i32, (*ptr.add(4))[0] as i32, (*ptr.add(3))[0] as i32, (*ptr.add(2))[0] as i32, (*ptr.add(1))[0] as i32, (*ptr.add(0))[0] as i32));
                    let v1 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[1] as i32, (*ptr.add(6))[1] as i32, (*ptr.add(5))[1] as i32, (*ptr.add(4))[1] as i32, (*ptr.add(3))[1] as i32, (*ptr.add(2))[1] as i32, (*ptr.add(1))[1] as i32, (*ptr.add(0))[1] as i32));
                    let v2 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[2] as i32, (*ptr.add(6))[2] as i32, (*ptr.add(5))[2] as i32, (*ptr.add(4))[2] as i32, (*ptr.add(3))[2] as i32, (*ptr.add(2))[2] as i32, (*ptr.add(1))[2] as i32, (*ptr.add(0))[2] as i32));
                    let v3 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[3] as i32, (*ptr.add(6))[3] as i32, (*ptr.add(5))[3] as i32, (*ptr.add(4))[3] as i32, (*ptr.add(3))[3] as i32, (*ptr.add(2))[3] as i32, (*ptr.add(1))[3] as i32, (*ptr.add(0))[3] as i32));

                    let (r0, r1, r2, r3) = extmul_8x(v0, v1, v2, v3, t_b0, t_b1, t_b2, t_b3);
                    let (l0, l1, l2, l3) = extmul_8x(v0, v1, v2, v3, omt_b0, omt_b1, omt_b2, omt_b3);

                    u642u32(eq_table.as_mut_ptr().add(i + current_len), r0, r1, r2, r3);
                    u642u32(eq_table.as_mut_ptr().add(i), l0, l1, l2, l3);
                }
            } else {
                let (left_half, right_half) = eq_table[..2 * current_len].split_at_mut(current_len);
                left_half.par_chunks_exact_mut(8).zip(right_half.par_chunks_exact_mut(8)).for_each(|(l_chunk, r_chunk)| {
                    let ptr = l_chunk.as_ptr();
                    let v0 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[0] as i32, (*ptr.add(6))[0] as i32, (*ptr.add(5))[0] as i32, (*ptr.add(4))[0] as i32, (*ptr.add(3))[0] as i32, (*ptr.add(2))[0] as i32, (*ptr.add(1))[0] as i32, (*ptr.add(0))[0] as i32));
                    let v1 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[1] as i32, (*ptr.add(6))[1] as i32, (*ptr.add(5))[1] as i32, (*ptr.add(4))[1] as i32, (*ptr.add(3))[1] as i32, (*ptr.add(2))[1] as i32, (*ptr.add(1))[1] as i32, (*ptr.add(0))[1] as i32));
                    let v2 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[2] as i32, (*ptr.add(6))[2] as i32, (*ptr.add(5))[2] as i32, (*ptr.add(4))[2] as i32, (*ptr.add(3))[2] as i32, (*ptr.add(2))[2] as i32, (*ptr.add(1))[2] as i32, (*ptr.add(0))[2] as i32));
                    let v3 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[3] as i32, (*ptr.add(6))[3] as i32, (*ptr.add(5))[3] as i32, (*ptr.add(4))[3] as i32, (*ptr.add(3))[3] as i32, (*ptr.add(2))[3] as i32, (*ptr.add(1))[3] as i32, (*ptr.add(0))[3] as i32));

                    let (r0, r1, r2, r3) = extmul_8x(v0, v1, v2, v3, t_b0, t_b1, t_b2, t_b3);
                    let (l0, l1, l2, l3) = extmul_8x(v0, v1, v2, v3, omt_b0, omt_b1, omt_b2, omt_b3);

                    u642u32(r_chunk.as_mut_ptr(), r0, r1, r2, r3);
                    u642u32(l_chunk.as_mut_ptr(), l0, l1, l2, l3);
                });
            }
        }
    }
    eq_table
}

pub unsafe fn build_eq_tbl_split(tau: &[[u32; 4]]) -> (Vec<[u32; 4]>, Vec<[u32; 4]>) {
    let high_tau = &tau[0..18]; 
    let low_tau = &tau[18..28];

    let eq_high = build_eq_tbl_base(high_tau);
    let eq_low = build_eq_tbl_base(low_tau);

    (eq_low, eq_high)
}

pub unsafe fn build_eq_tbl(tau: &[[u32; 4]]) -> Vec<[u32; 4]> {
    let n = tau.len();
    if n <= 14 {
        return build_eq_tbl_base(tau);
    }
    let k = 14;
    let low_tau = &tau[0..k];
    let high_tau = &tau[k..];
    let eq_low = build_eq_tbl_base(low_tau);
    let eq_high = build_eq_tbl_base(high_tau);
    let total_size = 1 << n;
    let mut eq_table: Vec<[u32; 4]> = Vec::with_capacity(total_size);
    eq_table.set_len(total_size);
    let chunk_size = 1 << k; // 16384
    eq_table.par_chunks_exact_mut(chunk_size).enumerate().for_each(|(i, chunk)| {
        let h_val = eq_high[i];
        let chunk_ptr = chunk.as_mut_ptr();
        let eq_low_ptr = eq_low.as_ptr();
        let h_b0 = _mm512_set1_epi64(h_val[0] as i64);
        let h_b1 = _mm512_set1_epi64(h_val[1] as i64);
        let h_b2 = _mm512_set1_epi64(h_val[2] as i64);
        let h_b3 = _mm512_set1_epi64(h_val[3] as i64);
        let mut j = 0;
        while j < chunk_size {
            let ptr = eq_low_ptr.add(j);
            let v0 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[0] as i32, (*ptr.add(6))[0] as i32, (*ptr.add(5))[0] as i32, (*ptr.add(4))[0] as i32, (*ptr.add(3))[0] as i32, (*ptr.add(2))[0] as i32, (*ptr.add(1))[0] as i32, (*ptr.add(0))[0] as i32));
            let v1 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[1] as i32, (*ptr.add(6))[1] as i32, (*ptr.add(5))[1] as i32, (*ptr.add(4))[1] as i32, (*ptr.add(3))[1] as i32, (*ptr.add(2))[1] as i32, (*ptr.add(1))[1] as i32, (*ptr.add(0))[1] as i32));
            let v2 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[2] as i32, (*ptr.add(6))[2] as i32, (*ptr.add(5))[2] as i32, (*ptr.add(4))[2] as i32, (*ptr.add(3))[2] as i32, (*ptr.add(2))[2] as i32, (*ptr.add(1))[2] as i32, (*ptr.add(0))[2] as i32));
            let v3 = _mm512_cvtepu32_epi64(_mm256_set_epi32((*ptr.add(7))[3] as i32, (*ptr.add(6))[3] as i32, (*ptr.add(5))[3] as i32, (*ptr.add(4))[3] as i32, (*ptr.add(3))[3] as i32, (*ptr.add(2))[3] as i32, (*ptr.add(1))[3] as i32, (*ptr.add(0))[3] as i32));
            let (r0, r1, r2, r3) = extmul_8x(v0, v1, v2, v3, h_b0, h_b1, h_b2, h_b3);
            u642u32(chunk_ptr.add(j), r0, r1, r2, r3);
            j += 8;
        }
    });

    eq_table
}