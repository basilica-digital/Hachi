use std::arch::x86_64::*;
use crate::math::field_simd::*;
use crate::Align64;
use crate::math::rns::*;
use tfhe_ntt::prime32::Plan;
use std::slice;
use crate::Q;

pub unsafe fn commit(u: &mut [u32], r: &mut [u32], t: &mut [u32], s: &[u8], acap: &[u32], bcap: &[u32]) { 
    let n = 1 << 10;
    let plans: Vec<Plan> = vec![
        Plan::try_new(1024, 759207937).unwrap(),
        Plan::try_new(1024, 759304193).unwrap(),
        Plan::try_new(2048, 2079301633).unwrap(),
        Plan::try_new(2048, 2079305729).unwrap(),
    ];
    let width = 1 << 10;
    let height = 1 << 13;
    let mut t_vec = vec![_mm512_set1_epi32(0); width * 64];
    let mut acap1_vec = vec![_mm512_set1_epi32(0); height * 64];
    let mut acap2_vec = vec![_mm512_set1_epi32(0); height * 64];
    let mut bcap1_vec = vec![_mm512_set1_epi32(0); height * 128];
    let mut bcap2_vec = vec![_mm512_set1_epi32(0); height * 128];
    let acap1_slice = std::slice::from_raw_parts_mut(acap1_vec.as_mut_ptr() as *mut u32, height * n);
    let acap2_slice = std::slice::from_raw_parts_mut(acap2_vec.as_mut_ptr() as *mut u32, height * n);
    let bcap1_slice = std::slice::from_raw_parts_mut(bcap1_vec.as_mut_ptr() as *mut u32, height*n*2);
    let bcap2_slice = std::slice::from_raw_parts_mut(bcap2_vec.as_mut_ptr() as *mut u32, height*n*2);
    // RNS and ntt on A, B
    for i in 0..height{
        for j in (0..1024).step_by(16) {
            rns_decompose_16_elements(&acap[(i<<10)+j..(i<<10)+(j+16)], &mut acap1_slice[(i<<10)+j..(i<<10)+(j+16)], &mut acap2_slice[(i<<10)+j..(i<<10)+(j+16)], 0);
            rns_2048(
                &bcap[(i<<10)+j .. (i<<10)+(j+16)], 
                &mut bcap1_slice[(i<<11)+j..(i<<11)+(j+16)], 
                &mut bcap2_slice[(i<<11)+j..(i<<11)+(j+16)]
            );
        }
        let mut a1 = &mut acap1_slice[i<<10..(i+1)<<10];
        let mut a2 = &mut acap2_slice[i<<10..(i+1)<<10];
        let mut b1 = &mut bcap1_slice[i<<11 .. (i+1)<<11];
        let mut b2 = &mut bcap2_slice[i<<11 .. (i+1)<<11];
        plans[0].fwd(&mut a1);
        plans[1].fwd(&mut a2);
        plans[2].fwd(&mut b1);
        plans[3].fwd(&mut b2);
    }
    let gt_ptr = t.as_mut_ptr() as *mut __m512i;
    let t_ptr = t_vec.as_mut_ptr();
    let u_ptr = u.as_mut_ptr() as *mut __m512i;
    let acap1_ptr = acap1_vec.as_mut_ptr();
    let acap2_ptr = acap2_vec.as_mut_ptr();
    let bcap1_ptr = bcap1_vec.as_mut_ptr();
    let bcap2_ptr = bcap2_vec.as_mut_ptr();
    // As
    let mut acc1 = vec![_mm512_set1_epi32(0); 1<<16];
    let mut acc2 = vec![_mm512_set1_epi32(0); 1<<16];
    // constants for montproduct
    let q1 = _mm512_set1_epi32(759207937);
    let q1_inv = _mm512_set1_epi32(754935809);
    let q2 = _mm512_set1_epi32(759304193);
    let q2_inv = _mm512_set1_epi32(331214849);
    let q3 = _mm512_set1_epi32(2079301633);
    let q3_inv = _mm512_set1_epi32(-1475321855);
    let q4 = _mm512_set1_epi32(2079305729);
    let q4_inv = _mm512_set1_epi32(-1659875327);
    let mask15 = _mm512_set1_epi32(15);
    let mask_0f = _mm_set1_epi8(0x0F);
    let acc1_ptr_base = acc1.as_mut_ptr() as *mut __m512i;
    let acc2_ptr_base = acc2.as_mut_ptr() as *mut __m512i;
    let acap1_ptr_simd = acap1_ptr as *const __m512i;
    let acap2_ptr_simd = acap2_ptr as *const __m512i;
    let s_ptr_packed = s.as_ptr() as *const u8;
    for j in 0..height {
        let s_row_offset_bytes = j << 19; 
        let a1_ptr = acap1_ptr_simd.add(j << 6);
        let a2_ptr = acap2_ptr_simd.add(j << 6);
        for i in 0..width{
            let mut s1 = Align64([0u32; 1024]);
            let mut s2 = Align64([0u32; 1024]);
            let s_raw_ptr = s_ptr_packed.add(s_row_offset_bytes + (i << 9)); 
            let s1_ptr = s1.0.as_mut_ptr() as *mut __m512i;
            let s2_ptr = s2.0.as_mut_ptr() as *mut __m512i;
            for k_step in 0..32 {
                let packed_128 = _mm_loadu_si128(s_raw_ptr.add(k_step * 16) as *const __m128i);
                let bytes_0 = _mm_and_si128(packed_128, mask_0f);
                let s_vec_0 = _mm512_cvtepu8_epi32(bytes_0);
                let bytes_1 = _mm_and_si128(_mm_srli_epi16(packed_128, 4), mask_0f);
                let s_vec_1 = _mm512_cvtepu8_epi32(bytes_1);
                _mm512_store_si512(s1_ptr.add(k_step * 2), s_vec_0);
                _mm512_store_si512(s2_ptr.add(k_step * 2), s_vec_0);
                _mm512_store_si512(s1_ptr.add(k_step * 2 + 1), s_vec_1);
                _mm512_store_si512(s2_ptr.add(k_step * 2 + 1), s_vec_1);
            }
            plans[0].fwd(&mut s1.0);
            plans[1].fwd(&mut s2.0);
            let s1_ptr = s1.0.as_ptr() as *const __m512i;
            let s2_ptr = s2.0.as_ptr() as *const __m512i;
            let acc1_ptr_i = acc1_ptr_base.add(i << 6);
            let acc2_ptr_i = acc2_ptr_base.add(i << 6);
            for k in 0..64 {
                let s1_vec = _mm512_load_si512(s1_ptr.add(k));
                let s2_vec = _mm512_load_si512(s2_ptr.add(k));
                let a1_vec = _mm512_load_si512(a1_ptr.add(k));
                let a2_vec = _mm512_load_si512(a2_ptr.add(k));
                let mut acc1_val = _mm512_load_si512(acc1_ptr_i.add(k));
                let mut acc2_val = _mm512_load_si512(acc2_ptr_i.add(k));
                let p1 = montproduct_759207937(s1_vec, a1_vec, q1, q1_inv);
                let p2 = montproduct_759304193(s2_vec, a2_vec, q2, q2_inv);
                
                let sum1 = _mm512_add_epi32(acc1_val, p1);
                let ge_q1 = _mm512_cmpge_epu32_mask(sum1, q1);
                acc1_val = _mm512_mask_sub_epi32(sum1, ge_q1, sum1, q1);
                
                let sum2 = _mm512_add_epi32(acc2_val, p2);
                let ge_q2 = _mm512_cmpge_epu32_mask(sum2, q2);
                acc2_val = _mm512_mask_sub_epi32(sum2, ge_q2, sum2, q2);
                
                _mm512_store_si512(acc1_ptr_i.add(k), acc1_val);
                _mm512_store_si512(acc2_ptr_i.add(k), acc2_val);
            }
        }
    }

    for i in 0..width{
        for j in 0..64{
            acc1[(i<<6)+j] = barrett_mul_4194304_759207937(acc1[(i<<6)+j]);
            acc2[(i<<6)+j] = barrett_mul_4194304_759304193(acc2[(i<<6)+j]);
        }
        let acc1_ptr = acc1.as_mut_ptr().add(i<<6) as *mut u32;
        let acc2_ptr = acc2.as_mut_ptr().add(i<<6) as *mut u32;
        plans[0].inv(slice::from_raw_parts_mut(acc1_ptr, 1024));
        plans[1].inv(slice::from_raw_parts_mut(acc2_ptr, 1024));
        let v1_ref = &*(acc1.as_ptr().add(i<<6) as *const [__m512i; 64]);
        let v2_ref = &*(acc2.as_ptr().add(i<<6) as *const [__m512i; 64]);
        irns(v1_ref, v2_ref, t_ptr, i<<6);
        // decompose t
        for j in 0..64{
            let mut val = _mm512_loadu_si512(t_ptr.add((i<<6)+j));
            for k in 0..8{
                _mm512_storeu_si512((gt_ptr.add((i<<9)+(k<<6)+j)) as *mut _ , _mm512_and_epi32(mask15, val));
                val = _mm512_srli_epi32(val, 4);
            }
        }
    }

    // Bt
    let mut v1 = [_mm512_set1_epi32(0); 128];
    let mut v2 = [_mm512_set1_epi32(0); 128];
    // Ring inner product
    for i in 0..height{
        let mut gt1 = Align64([0u32; 2048]);
        gt1.0[0..1024].copy_from_slice(&t[i<<10..(i+1)<<10]);
        let mut gt2 = gt1;
        plans[2].fwd(&mut gt1.0);
        plans[3].fwd(&mut gt2.0);
        let gt1_ptr = gt1.0.as_mut_ptr() as *mut __m512i;
        let gt2_ptr = gt2.0.as_mut_ptr() as *mut __m512i;
        for j in 0..128 {
            let gt1_vec = _mm512_loadu_si512(gt1_ptr.add(j));
            let gt2_vec = _mm512_loadu_si512(gt2_ptr.add(j));
            let b1_vec = _mm512_loadu_si512(bcap1_ptr.add((i<<7)+j));
            let b2_vec = _mm512_loadu_si512(bcap2_ptr.add((i<<7)+j));
            v1[j] = barrett_fake_2079301633(_mm512_add_epi32(v1[j], montproduct(gt1_vec, b1_vec, q3, q3_inv)));
            v2[j] = barrett_fake_2079305729(_mm512_add_epi32(v2[j], montproduct(gt2_vec, b2_vec, q4, q4_inv)));
        }
    }
    for i in 0..128 {
        v1[i] = barrett_mul_2097152_2079301633(v1[i]);
        v2[i] = barrett_mul_2097152_2079305729(v2[i]);
    }
    plans[2].inv(slice::from_raw_parts_mut(v1.as_mut_ptr() as *mut u32, 2048));
    plans[3].inv(slice::from_raw_parts_mut(v2.as_mut_ptr() as *mut u32, 2048));
    let mut full_res = Align64([0u32; 2048]);
    irns_2048(&v1, &v2, full_res.0.as_mut_ptr() as *mut __m512i, 0);
    let q_vec = _mm512_set1_epi32(((1u64<<32)-99) as i32);
    let full_res_ptr = full_res.0.as_ptr() as *const __m512i;
    let r_ptr = r.as_mut_ptr() as *mut __m512i;
    let u_ptr_simd = u.as_mut_ptr() as *mut __m512i;
    for j in 0..64 {
        let low_vec = _mm512_load_si512(full_res_ptr.add(j));
        let high_vec = _mm512_load_si512(full_res_ptr.add(64 + j)); 
        _mm512_storeu_si512(r_ptr.add(j), high_vec);
        let lt_mask = _mm512_cmplt_epu32_mask(low_vec, high_vec);
        let diff_vec = _mm512_sub_epi32(low_vec, high_vec);
        let u_vec = _mm512_mask_add_epi32(diff_vec, lt_mask, diff_vec, q_vec);
        _mm512_storeu_si512(u_ptr_simd.add(j), u_vec);
    }
}