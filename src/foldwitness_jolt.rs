#[inline(always)]
pub unsafe fn jolt_mod_q(out: &mut [i32; 8], x: __m512i) {
    let mask32 = _mm512_set1_epi64(0xFFFFFFFF_i64);
    let mask31 = _mm512_set1_epi64(0x7FFFFFFF_i64);
    let c38    = _mm512_set1_epi64(38_i64);
    let c19    = _mm512_set1_epi64(19_i64);
    let q_vec  = _mm512_set1_epi64(0x7FFFFFED_i64);

    let a = _mm512_srli_epi64(x, 32);
    let b = _mm512_and_si512(x, mask32);
    let a38 = _mm512_mul_epu32(a, c38);
    let y = _mm512_add_epi64(a38, b);

    let c = _mm512_srli_epi64(y, 32);
    let d = _mm512_and_si512(y, mask32);
    let c38_2 = _mm512_mul_epu32(c, c38);
    let z = _mm512_add_epi64(c38_2, d);

    let e = _mm512_srli_epi64(z, 31);
    let f = _mm512_and_si512(z, mask31);
    let e19 = _mm512_mul_epu32(e, c19);
    let w = _mm512_add_epi64(e19, f);

    let cmp = _mm512_cmp_epu64_mask(w, q_vec, 5);
    let result_64 = _mm512_mask_sub_epi64(w, cmp, w, q_vec);
    let result_32 = _mm512_cvtepi64_epi32(result_64);
    _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, result_32);
}

pub unsafe fn sparse_poly_index_jolt(id:&mut [u8], c: &[u64]){
    for i in 0..4{
        for j in 0..4{
            let val = c[i*4+j];
            if val != 0 {
                id[i] = (j*64 + val.trailing_zeros() as usize) as u8;
                break;
            }
        }
    }
}

#[inline(always)]
pub unsafe fn fold_witness_row_jolt_4(z: &mut [i32], id: &[u8], c: &[i32]) {
    const batch: usize = 4;
    // batch_stride = 257*8;
    
    let mut res = vec![_mm512_setzero_si512(); 8224];
    let res_ptr = res.as_mut_ptr();
    let c_ptr = c.as_ptr() as *const i32;
    let id_ptr = id.as_ptr();

    let mut lookuptable = [0u16; 256];
    for j in 0..256 {
        lookuptable[j] = ((j&7)*257 + (j>>3)) as u16;
    }

    for i in 0..(1 << 16) {
        let r = i*4;
        let row_c_ptr = c_ptr.add(i*1024);
        let mut p_stack = [std::ptr::null_mut::<__m512i>(); 16];

        for l in 0..4 {
            let ids_u32 = *(id_ptr.add((l << 18) + r) as *const u32);
            let b_off = l*2056;

            let j0 = (ids_u32 & 0xFF) as usize;
            let j1 = ((ids_u32 >> 8) & 0xFF) as usize;
            let j2 = ((ids_u32 >> 16) & 0xFF) as usize;
            let j3 = (ids_u32 >> 24) as usize;

            p_stack[l*4+0] = res_ptr.add(b_off + lookuptable[j0] as usize);
            p_stack[l*4+1] = res_ptr.add(b_off + lookuptable[j1] as usize + 32);
            p_stack[l*4+2] = res_ptr.add(b_off + lookuptable[j2] as usize + 64);
            p_stack[l*4+3] = res_ptr.add(b_off + lookuptable[j3] as usize + 96);
        }

        let p00 = p_stack[0];
        let p01 = p_stack[1];
        let p02 = p_stack[2];
        let p03 = p_stack[3];
        let p04 = p_stack[4];
        let p05 = p_stack[5];
        let p06 = p_stack[6];
        let p07 = p_stack[7];
        let p08 = p_stack[8];
        let p09 = p_stack[9];
        let p010 = p_stack[10];
        let p011 = p_stack[11];
        let p012 = p_stack[12];
        let p013 = p_stack[13];
        let p014 = p_stack[14];
        let p015 = p_stack[15];

        for k in (0..128).step_by(2) {
            let v512 = _mm512_loadu_si512(row_c_ptr.add(k * 8) as *const __m512i);
            let cv0 = _mm512_cvtepi32_epi64(_mm512_castsi512_si256(v512));
            let cv1 = _mm512_cvtepi32_epi64(_mm512_extracti64x4_epi64::<1>(v512));
            
            // 1
            let p0 = p00.add(k);
            let p1 = p01.add(k);
            let p2 = p02.add(k);
            let p3 = p03.add(k);
            *p0 = _mm512_add_epi64(cv0, *p0);
            *(p0.add(1)) = _mm512_add_epi64(cv1, *(p0.add(k+1)));
            *p1 = _mm512_add_epi64(cv0, *p1);
            *(p1.add(1)) = _mm512_add_epi64(cv1, *(p1.add(1)));
            *p2 = _mm512_add_epi64(cv0, *p2);
            *(p2.add(1)) = _mm512_add_epi64(cv1, *(p2.add(1)));
            *p3 = _mm512_add_epi64(cv0, *p3);
            *(p3.add(1)) = _mm512_add_epi64(cv1, *(p3.add(1)));
            
            // 2
            let p0 = p04.add(k);
            let p1 = p05.add(k);
            let p2 = p06.add(k);
            let p3 = p07.add(k);
            *p0 = _mm512_add_epi64(cv0, *p0);
            *(p0.add(1)) = _mm512_add_epi64(cv1, *(p0.add(1)));
            *p1 = _mm512_add_epi64(cv0, *p1);
            *(p1.add(1)) = _mm512_add_epi64(cv1, *(p1.add(1)));
            *p2 = _mm512_add_epi64(cv0, *p2);
            *(p2.add(1)) = _mm512_add_epi64(cv1, *(p2.add(1)));
            *p3 = _mm512_add_epi64(cv0, *p3);
            *(p3.add(1)) = _mm512_add_epi64(cv1, *(p3.add(1)));
            
            // 3
            let p0 = p08.add(k);
            let p1 = p09.add(k);
            let p2 = p010.add(k);
            let p3 = p011.add(k);
            *p0 = _mm512_add_epi64(cv0,*p0);
            *(p0.add(1)) = _mm512_add_epi64(cv1, *(p0.add(1)));
            *p1 = _mm512_add_epi64(cv0, *p1);
            *(p1.add(1)) = _mm512_add_epi64(cv1, *(p1.add(1)));
            *p2 = _mm512_add_epi64(cv0, *p2);
            *(p2.add(1)) = _mm512_add_epi64(cv1, *(p2.add(1)));
            *p3 = _mm512_add_epi64(cv0, *p3);
            *(p3.add(1)) = _mm512_add_epi64(cv1, *(p3.add(1)));
            
            // 4
            let p0 = p012.add(k);
            let p1 = p013.add(k);
            let p2 = p014.add(k);
            let p3 = p015.add(k);
            *p0 = _mm512_add_epi64(cv0, *p0);
            *(p0.add(1)) = _mm512_add_epi64(cv1, *(p0.add(1)));
            *p1 = _mm512_add_epi64(cv0, *p1);
            *(p1.add(1)) = _mm512_add_epi64(cv1, *(p1.add(1)));
            *p2 = _mm512_add_epi64(cv0, *p2);
            *(p2.add(1)) = _mm512_add_epi64(cv1, *(p2.add(1)));
            *p3 = _mm512_add_epi64(cv0, *p3);
            *(p3.add(1)) = _mm512_add_epi64(cv1, *(p3.add(1)));
        
        }
    }
    const BIN_SIZE: usize = 257;
    // merge
    for i in 0..batch {
        let b_off = i*2056;
        let mut p1 = _mm512_setzero_si512(); let mut p2 = _mm512_setzero_si512();
        let mut p3 = _mm512_setzero_si512(); let mut p4 = _mm512_setzero_si512();
        let mut p5 = _mm512_setzero_si512(); let mut p6 = _mm512_setzero_si512();
        let mut p7 = _mm512_setzero_si512();

        for k in 0..128 {
            let mut acc = *res_ptr.add(b_off + 0 * BIN_SIZE + k);
            let c1 = *res_ptr.add(b_off + 1 * BIN_SIZE + k);
            let c2 = *res_ptr.add(b_off + 2 * BIN_SIZE + k);
            let c3 = *res_ptr.add(b_off + 3 * BIN_SIZE + k);
            let c4 = *res_ptr.add(b_off + 4 * BIN_SIZE + k);
            let c5 = *res_ptr.add(b_off + 5 * BIN_SIZE + k);
            let c6 = *res_ptr.add(b_off + 6 * BIN_SIZE + k);
            let c7 = *res_ptr.add(b_off + 7 * BIN_SIZE + k);

            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<7>(c1, p1)); p1 = c1;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<6>(c2, p2)); p2 = c2;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<5>(c3, p3)); p3 = c3;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<4>(c4, p4)); p4 = c4;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<3>(c5, p5)); p5 = c5;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<2>(c6, p6)); p6 = c6;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<1>(c7, p7)); p7 = c7;
            *res_ptr.add(b_off + 0 * BIN_SIZE + k) = acc;
        }

        for k in 128..256 {
            let mut acc = *res_ptr.add(b_off + 0 * BIN_SIZE + k);
            let c1 = *res_ptr.add(b_off + 1 * BIN_SIZE + k);
            let c2 = *res_ptr.add(b_off + 2 * BIN_SIZE + k);
            let c3 = *res_ptr.add(b_off + 3 * BIN_SIZE + k);
            let c4 = *res_ptr.add(b_off + 4 * BIN_SIZE + k);
            let c5 = *res_ptr.add(b_off + 5 * BIN_SIZE + k);
            let c6 = *res_ptr.add(b_off + 6 * BIN_SIZE + k);
            let c7 = *res_ptr.add(b_off + 7 * BIN_SIZE + k);

            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<7>(c1, p1)); p1 = c1;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<6>(c2, p2)); p2 = c2;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<5>(c3, p3)); p3 = c3;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<4>(c4, p4)); p4 = c4;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<3>(c5, p5)); p5 = c5;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<2>(c6, p6)); p6 = c6;
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<1>(c7, p7)); p7 = c7;
            let target_ptr = res_ptr.add(b_off + 0 * BIN_SIZE + (k - 128));
            *target_ptr = _mm512_sub_epi64(*target_ptr, acc);
        }

        let offset_q = _mm512_set1_epi64(((1i64 << 31) - 19) << 19);
        let z_base_ptr = z.as_mut_ptr().add(i * 1024);
        for k in 0..128 {
            let val = _mm512_add_epi64(*res_ptr.add(b_off + 0 * BIN_SIZE + k), offset_q);
            jolt_mod_q(&mut *(z_base_ptr.add(k * 8) as *mut [i32; 8]), val);
        }
    }
}

#[inline(never)]
pub unsafe fn fold_witness_jolt_4(z: &mut [i32], c: &[i32], s: &[u8]) {
    z.par_chunks_mut(1 << 12).enumerate().for_each(|(i, z_row)| {
        fold_witness_row_jolt_4(z_row, &s[i*((1<<20))..], c);
    });
}