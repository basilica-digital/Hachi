pub unsafe fn sparse_poly_index(id:&mut [i16], c: &[i16]){
    let mut count = 0;
    for i in 0..1024{
        if c[i] == 1{
            id[count] = (i+1) as i16;
            count = count+1;
        }else if c[i] == -1{
            id[count] = -((i+1) as i16);
            count = count+1;
        }
        if count == 16{
            break;
        }
    }
}

pub unsafe fn sparse_poly_index_u32(id:&mut [i16], c: &[u32]){
    let q = (1u64 << 32) - 99;
    let mut count = 0;
    for i in 0..1024{
        if c[i] == 1{
            id[count] = (i+1) as i16;
            count = count+1;
        }else if c[i] == (q-1) as u32 {
            id[count] = -((i+1) as i16);
            count = count+1;
        }
        if count == 16{
            break;
        }
    }
}

#[inline(always)]
pub unsafe fn rotate_and_add(s: &[i16], i: usize, z: &mut [i16]) {
    let n = 1024;
    let src1 = &s[n-i..n];
    let dst1 = &mut z[0..i];
    vector_sub_assign(src1, dst1);

    let src2 = &s[0..n - i];
    let dst2 = &mut z[i..n];
    vector_add_assign(src2, dst2);
}

#[inline(always)]
pub unsafe fn rotate_and_add_u32(s: &[u32], i: usize, z: &mut [u32]) {
    let n = 1024;
    let src1 = &s[n-i..n];
    let dst1 = &mut z[0..i];
    vector_sub_assign_u32(src1, dst1);

    let src2 = &s[0..n - i];
    let dst2 = &mut z[i..n];
    vector_add_assign_u32(src2, dst2);
}

#[inline(always)]
pub unsafe fn rotate_and_sub(s: &[i16], i: usize, z: &mut [i16]) {
    let n = 1024;
    let src1 = &s[n-i..n];
    let dst1 = &mut z[0..i];
    vector_add_assign(src1, dst1);

    let src2 = &s[0..n-i];
    let dst2 = &mut z[i..n];
    vector_sub_assign(src2, dst2);
}

#[inline(always)]
pub unsafe fn rotate_and_sub_u32(s: &[u32], i: usize, z: &mut [u32]) {
    let n = 1024;
    let src1 = &s[n-i..n];
    let dst1 = &mut z[0..i];
    vector_add_assign_u32(src1, dst1);

    let src2 = &s[0..n-i];
    let dst2 = &mut z[i..n];
    vector_sub_assign_u32(src2, dst2);
}

#[inline(always)]
unsafe fn vector_add_assign(src: &[i16], dst: &mut [i16]) {
    let mut len = src.len();
    let mut s_ptr = src.as_ptr();
    let mut d_ptr = dst.as_mut_ptr();

    while len >= 32 {
        let v_s = _mm512_loadu_si512(s_ptr as *const _);
        let v_d = _mm512_loadu_si512(d_ptr as *const _);
        let res = _mm512_add_epi16(v_d, v_s);
        _mm512_storeu_si512(d_ptr as *mut _, res);
        s_ptr = s_ptr.add(32); d_ptr = d_ptr.add(32); len -= 32;
    }
    if len > 0 {
        let mask = (1u32 << len) - 1;
        let k = _cvtu32_mask32(mask);
        let v_s = _mm512_maskz_loadu_epi16(k, s_ptr as *const _);
        let v_d = _mm512_maskz_loadu_epi16(k, d_ptr as *const _);
        let res = _mm512_add_epi16(v_d, v_s);
        _mm512_mask_storeu_epi16(d_ptr as *mut _, k, res);
    }
}

#[inline(always)]
unsafe fn vector_add_assign_u32(src: &[u32], dst: &mut [u32]) {
    let mut len = src.len();
    let mut s_ptr = src.as_ptr();
    let mut d_ptr = dst.as_mut_ptr();
    let q = (1u64 << 32) - 99;
    let q_vec = _mm512_set1_epi32(q as i32);

    while len >= 16 {
        let v_s = _mm512_loadu_si512(s_ptr as *const _);
        let v_d = _mm512_loadu_si512(d_ptr as *const _);
        let mut res = _mm512_add_epi32(v_d, v_s);
        
        // Modular reduction: if res >= q, res -= q
        let mask = _mm512_cmpge_epu32_mask(res, q_vec);
        res = _mm512_mask_sub_epi32(res, mask, res, q_vec);
        
        _mm512_storeu_si512(d_ptr as *mut _, res);
        s_ptr = s_ptr.add(16); d_ptr = d_ptr.add(16); len -= 16;
    }
    for i in 0..len {
        dst[i] = ((dst[i] as u64 + src[i] as u64) % q) as u32;
    }
}

#[inline(always)]
unsafe fn vector_sub_assign(src: &[i16], dst: &mut [i16]) {
    let mut len = src.len();
    let mut s_ptr = src.as_ptr();
    let mut d_ptr = dst.as_mut_ptr();

    while len >= 32 {
        let v_s = _mm512_loadu_si512(s_ptr as *const _);
        let v_d = _mm512_loadu_si512(d_ptr as *const _);
        let res = _mm512_sub_epi16(v_d, v_s);
        _mm512_storeu_si512(d_ptr as *mut _, res);
        s_ptr = s_ptr.add(32); d_ptr = d_ptr.add(32); len -= 32;
    }
    if len > 0 {
        let mask = (1u32 << len) - 1;
        let k = _cvtu32_mask32(mask);
        let v_s = _mm512_maskz_loadu_epi16(k, s_ptr as *const _);
        let v_d = _mm512_maskz_loadu_epi16(k, d_ptr as *const _);
        let res = _mm512_sub_epi16(v_d, v_s);
        _mm512_mask_storeu_epi16(d_ptr as *mut _, k, res);
    }
}

#[inline(always)]
unsafe fn vector_sub_assign_u32(src: &[u32], dst: &mut [u32]) {
    let mut len = src.len();
    let mut s_ptr = src.as_ptr();
    let mut d_ptr = dst.as_mut_ptr();
    let q = (1u64 << 32) - 99;
    let q_vec = _mm512_set1_epi32(q as i32);

    while len >= 16 {
        let v_s = _mm512_loadu_si512(s_ptr as *const _);
        let v_d = _mm512_loadu_si512(d_ptr as *const _);
        
        // dst = (dst - src + q) % q
        // Use mask to handle underflow
        let mask = _mm512_cmplt_epu32_mask(v_d, v_s);
        let mut res = _mm512_sub_epi32(v_d, v_s);
        res = _mm512_mask_add_epi32(res, mask, res, q_vec);
        
        _mm512_storeu_si512(d_ptr as *mut _, res);
        s_ptr = s_ptr.add(16); d_ptr = d_ptr.add(16); len -= 16;
    }
    for i in 0..len {
        dst[i] = ((dst[i] as u64 + q - src[i] as u64) % q) as u32;
    }
}

#[inline(always)]
pub unsafe fn fold_witness_ring(z: &mut [i16], id: &[i16], s: &[i16]){
    for i in 0..16 {
        if id[i] > 0 {
            let index = (id[i] - 1) as usize;
            rotate_and_add(s, index, z);
        } else {
            let index = ((-id[i]) - 1) as usize;
            rotate_and_sub(s, index, z);
        }
    }
}

#[inline(always)]
pub unsafe fn fold_witness_ring_u32(z: &mut [u32], id: &[i16], s: &[u32]){
    for i in 0..16 {
        if id[i] > 0 {
            let index = (id[i] - 1) as usize;
            rotate_and_add_u32(s, index, z);
        } else if id[i] < 0 {
            let index = ((-id[i]) - 1) as usize;
            rotate_and_sub_u32(s, index, z);
        }
    }
}

// #[inline(always)]
// pub unsafe fn fold_witness_ring_test(z: &mut [i16], id: &[i16], s: &[i16]){
//     let mut large: [i16;2048] = [0i16;2048];
//     for i in 0..16 {
//         if id[i] > 0 {
//             let index = (id[i] - 1) as usize;
//             for j in 0..1024{
//                 large[index+j] = large[index+j] + s[j];
//             }
//         } else {
//             let index = ((-id[i]) - 1) as usize;
//             for j in 0..1024{
//                 large[index+j] = large[index+j] - s[j];
//             }
//         }
//     }
//     for i in 0..1024{
//         z[i] = z[i] + large[i] - large[i+1024];
//     }
// }

#[inline(always)]
pub unsafe fn fold_witness_row(z: &mut [i16], id: &[i16], s: &[i16]){
    for i in 0..1024{
        fold_witness_ring(z, &id[i*16..], &s[i*(1<<10)..]);
        //fold_witness_ring_test(z, &id[i*16..], &s[i*(1<<10)..]);
    }
}

#[inline(always)]
pub unsafe fn fold_witness_row_u32(z: &mut [u32], id: &[i16], s: &[u32]){
    for i in 0..1024{
        fold_witness_ring_u32(z, &id[i*16..], &s[i*(1<<10)..]);
    }
}

#[inline(never)]
pub unsafe fn fold_witness(z: &mut [i16], c: &[i16], s: &[i16]) {
    let mut id = vec![0i16; 16384];
    id.par_chunks_mut(16).enumerate().for_each(|(i, id_slice)| {
        sparse_poly_index(id_slice, &c[1024 * i..]);
    });
    z.par_chunks_mut(1 << 10).enumerate().for_each(|(i, z_row)| {
        fold_witness_row(z_row, &id, &s[i*(1 << 20)..]);
    });
}
#[inline(never)]
pub unsafe fn fold_witness_u32(z: &mut [u32], c: &[u32], s: &[u32]) {
    let mut id = vec![0i16; 16384];
    id.par_chunks_mut(16).enumerate().for_each(|(i, id_slice)| {
        sparse_poly_index_u32(id_slice, &c[1024 * i..]);
    });
    z.par_chunks_mut(1 << 10).enumerate().for_each(|(i, z_row)| {
        fold_witness_row_u32(z_row, &id, &s[i*(1 << 20)..]);
    });
}
