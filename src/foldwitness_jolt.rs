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
    *out = std::mem::transmute::<__m256i, [i32; 8]>(result_32);
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
pub unsafe fn fold_witness_row_jolt(z: &mut [i32], id: &[u8], c: &[i32], batch: usize){
    let mut res: Vec<[[__m512i; 256]; 8]> = vec![[[_mm512_setzero_si512(); 256]; 8]; batch];
    let c_ptr = c.as_ptr() as *const __m256i;
    // inner product
    for i in 0..(1<<16){
        // id
        for j in 0..4{
            // c 
            for k in 0..128{
                let v256 = _mm256_load_si256(c_ptr.add(i*128 + k));
                let c_vec = _mm512_cvtepi32_epi64(v256);
                // batch 
                for l in 0..batch{
                    let offset = id[l*(1<<18)+i*4+j] as usize;
                    let channel = offset % 8;
                    let base_idx = (offset / 8) + j*32;
                    res[l][channel][base_idx + k] = _mm512_add_epi64(c_vec, res[l][channel][base_idx + k]);
                }
            }
        }
    }
    // merge res
    for i in 0..batch{
        let mut prev_1 = _mm512_setzero_si512();
        let mut prev_2 = _mm512_setzero_si512();
        let mut prev_3 = _mm512_setzero_si512();
        let mut prev_4 = _mm512_setzero_si512();
        let mut prev_5 = _mm512_setzero_si512();
        let mut prev_6 = _mm512_setzero_si512();
        let mut prev_7 = _mm512_setzero_si512();

        for k in 0..128 {
            let mut acc = res[i][0][k];
            
            // j = 1, shift left by 1 element (64 bits), extract starting from 8-1 = 7
            let c1 = res[i][1][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<7>(c1, prev_1));
            prev_1 = c1;

            // j = 2, extract starting from 8-2 = 6
            let c2 = res[i][2][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<6>(c2, prev_2));
            prev_2 = c2;

            // j = 3, extract starting from 8-3 = 5
            let c3 = res[i][3][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<5>(c3, prev_3));
            prev_3 = c3;

            // j = 4, extract starting from 8-4 = 4
            let c4 = res[i][4][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<4>(c4, prev_4));
            prev_4 = c4;

            // j = 5, extract starting from 8-5 = 3
            let c5 = res[i][5][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<3>(c5, prev_5));
            prev_5 = c5;

            // j = 6, extract starting from 8-6 = 2
            let c6 = res[i][6][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<2>(c6, prev_6));
            prev_6 = c6;

            // j = 7, extract starting from 8-7 = 1
            let c7 = res[i][7][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<1>(c7, prev_7));
            prev_7 = c7;
            res[i][0][k] = acc;
        }

        for k in 128..256 {
            let mut acc = res[i][0][k];
            
            // j = 1, shift left by 1 element (64 bits), extract starting from 8-1 = 7
            let c1 = res[i][1][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<7>(c1, prev_1));
            prev_1 = c1;

            // j = 2, extract starting from 8-2 = 6
            let c2 = res[i][2][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<6>(c2, prev_2));
            prev_2 = c2;

            // j = 3, extract starting from 8-3 = 5
            let c3 = res[i][3][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<5>(c3, prev_3));
            prev_3 = c3;

            // j = 4, extract starting from 8-4 = 4
            let c4 = res[i][4][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<4>(c4, prev_4));
            prev_4 = c4;

            // j = 5, extract starting from 8-5 = 3
            let c5 = res[i][5][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<3>(c5, prev_5));
            prev_5 = c5;

            // j = 6, extract starting from 8-6 = 2
            let c6 = res[i][6][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<2>(c6, prev_6));
            prev_6 = c6;

            // j = 7, extract starting from 8-7 = 1
            let c7 = res[i][7][k];
            acc = _mm512_add_epi64(acc, _mm512_alignr_epi64::<1>(c7, prev_7));
            prev_7 = c7;
            res[i][0][k-128] = _mm512_sub_epi64(res[i][0][k-128], acc);
        }
    }
    for j in 0..batch {
        for i in 0..128{
            jolt_mod_q((&mut z[j*1024+i*8..j*1024+i*8+8]).try_into().unwrap(), res[j][0][i]);
        }
    }
}

#[inline(never)]
pub unsafe fn fold_witness_jolt(z: &mut [i32], c: &[i32], s: &[u64], batch: usize) {
    // id size: 2 + 16 + 6
    let mut id = vec![0u8; (1<<24)];
    id.par_chunks_mut(4).enumerate().for_each(|(i, id_slice)| {
        sparse_poly_index_jolt(id_slice, &s[16 * i..]);
    });
    z.par_chunks_mut((1 << 10)*batch).enumerate().for_each(|(i, z_row)| {
        fold_witness_row_jolt(z_row, &id[i*((1<<18)*batch)..], c, batch);
    });
}