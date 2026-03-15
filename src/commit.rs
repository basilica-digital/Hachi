#[inline(never)]
pub unsafe fn commit(u: &mut [u32], s: &mut [u32], acap: &mut [u32], bcap: &mut [u32]){
    let n = 1 << 10;
    let moduli = [759207937, 759304193];
    let plans: Vec<Plan> = moduli.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();
    let width = 1 << 10;
    let height = 1 << 13;

    let mut t_vec = vec![_mm512_set1_epi32(0); width * 64];
    let mut gt_vec = vec![_mm512_set1_epi32(0); height * 64];
    let mut acap1_vec = vec![_mm512_set1_epi32(0); height * 64];
    let mut acap2_vec = vec![_mm512_set1_epi32(0); height * 64];
    let mut bcap1_vec = vec![_mm512_set1_epi32(0); height * 64];
    let mut bcap2_vec = vec![_mm512_set1_epi32(0); height * 64];
    let acap1_slice = std::slice::from_raw_parts_mut(acap1_vec.as_mut_ptr() as *mut u32, height * n);
    let acap2_slice = std::slice::from_raw_parts_mut(acap2_vec.as_mut_ptr() as *mut u32, height * n);
    let bcap1_slice = std::slice::from_raw_parts_mut(bcap1_vec.as_mut_ptr() as *mut u32, height * n);
    let bcap2_slice = std::slice::from_raw_parts_mut(bcap2_vec.as_mut_ptr() as *mut u32, height * n);
    
    // RNS and ntt on A, B
    for i in 0..height{
        for j in (0..1024).step_by(16) {
            rns_decompose_16_elements(&acap[(i<<10)+j..(i<<10)+(j+16)], &mut acap1_slice[(i<<10)+j..(i<<10)+(j+16)], &mut acap2_slice[(i<<10)+j..(i<<10)+(j+16)], j);
            rns_decompose_16_elements(&bcap[(i<<10)+j..(i<<10)+(j+16)], &mut bcap1_slice[(i<<10)+j..(i<<10)+(j+16)], &mut bcap2_slice[(i<<10)+j..(i<<10)+(j+16)], j);
        }
        let mut a1 = &mut acap1_slice[i<<10..(i+1)<<10];
        let mut a2 = &mut acap2_slice[i<<10..(i+1)<<10];
        let mut b1 = &mut bcap1_slice[i<<10..(i+1)<<10];
        let mut b2 = &mut bcap2_slice[i<<10..(i+1)<<10];
        plans[0].fwd(&mut a1);
        plans[1].fwd(&mut a2);
        plans[0].fwd(&mut b1);
        plans[1].fwd(&mut b2);
    }

    let gt_ptr = gt_vec.as_mut_ptr();
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
    let mask15 = _mm512_set1_epi32(15);

    let mut s1 = Align64([0u32; 1024]);
    let mut s2 = Align64([0u32; 1024]);

    let acc1_ptr_base = acc1.as_mut_ptr() as *mut __m512i;
    let acc2_ptr_base = acc2.as_mut_ptr() as *mut __m512i;
    let acap1_ptr_simd = acap1_ptr as *const __m512i;
    let acap2_ptr_simd = acap2_ptr as *const __m512i;

    for j in 0..height {
        let original_j = j/8;
        let bit_shift = (j%8)*4;
        let shift_count = _mm_cvtsi32_si128(bit_shift as i32);
        let s_row_offset = original_j << 20;
        let a1_ptr = acap1_ptr_simd.add(j << 6);
        let a2_ptr = acap2_ptr_simd.add(j << 6);
        for i in 0..width{
            let mut s1 = Align64([0u32; 1024]);
            let mut s2 = Align64([0u32; 1024]);
            let s_raw_slice = &s[(s_row_offset + (i << 10)) ..]; 
            let s_raw_ptr = s_raw_slice.as_ptr() as *const __m512i;
            let s1_ptr = s1.0.as_mut_ptr() as *mut __m512i;
            let s2_ptr = s2.0.as_mut_ptr() as *mut __m512i;
            for k in 0..64 {
                let val = _mm512_loadu_si512(s_raw_ptr.add(k));
                let shifted = _mm512_srl_epi32(val, shift_count);
                let decomp = _mm512_and_epi32(shifted, mask15);
                _mm512_store_si512(s1_ptr.add(k), decomp);
                _mm512_store_si512(s2_ptr.add(k), decomp);
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
        irns(v1_ref, v2_ref, t_ptr, i);
    }

    // decompose t
    for i in 0..width{
        for j in 0..64{
            let mut val = _mm512_loadu_si512(t_ptr.add((i<<6)+j));
            for k in 0..8{
                _mm512_storeu_si512((gt_ptr.add((i<<9)+(k<<6)+j)) as *mut _ , _mm512_and_epi32(mask15, val));
                val = _mm512_srli_epi32(val, 4);
            }
        }
    }

    // Bt
    let mut v1 = [_mm512_set1_epi32(0); 64];
    let mut v2 = [_mm512_set1_epi32(0); 64];
    let gt_slice = std::slice::from_raw_parts(gt_ptr as *const u32, height * n);
    // Ring wise inner product
    for i in 0..height{
        let mut gt1 = Align64([0u32; 1024]);
        gt1.0.copy_from_slice(&gt_slice[i<<10..(i+1)<<10]);
        let mut gt2 = gt1;
        plans[0].fwd(&mut gt1.0);
        plans[1].fwd(&mut gt2.0);
        let gt1_ptr = gt1.0.as_mut_ptr() as *mut __m512i;
        let gt2_ptr = gt2.0.as_mut_ptr() as *mut __m512i;
        for j in 0..64{
            let gt1_vec = _mm512_loadu_si512(gt1_ptr.add(j));
            let gt2_vec = _mm512_loadu_si512(gt2_ptr.add(j));
            let b1_vec = _mm512_loadu_si512(bcap1_ptr.add((i<<6)+j));
            let b2_vec = _mm512_loadu_si512(bcap2_ptr.add((i<<6)+j));
            v1[j] = barrett_fake_759207937(_mm512_add_epi32(v1[j], montproduct_759207937(gt1_vec, b1_vec, q1, q1_inv)));
            v2[j] = barrett_fake_759304193(_mm512_add_epi32(v2[j], montproduct_759304193(gt2_vec, b2_vec, q2, q2_inv)));  
        }
    }
    for i in 0..64{
        v1[i] = barrett_mul_4194304_759207937(v1[i]);
        v2[i] = barrett_mul_4194304_759304193(v2[i]);
    }
    plans[0].inv(slice::from_raw_parts_mut(v1.as_mut_ptr() as *mut u32, 1024));
    plans[1].inv(slice::from_raw_parts_mut(v2.as_mut_ptr() as *mut u32, 1024));
    irns(&v1, &v2, u_ptr, 0);
}