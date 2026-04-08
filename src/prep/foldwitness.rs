use std::ptr;
use std::arch::x86_64::*;
use rayon::prelude::*;

#[inline(always)]
unsafe fn load_unpacked_s(s_ptr: *const u8) -> __m512i {
    let mask_0f = _mm_set1_epi8(0x0F);
    let packed_128 = _mm_loadu_si128(s_ptr as *const __m128i);
    let bytes_0 = _mm_and_si128(packed_128, mask_0f);
    let s_lo = _mm256_cvtepu8_epi16(bytes_0);
    let bytes_1 = _mm_and_si128(_mm_srli_epi16(packed_128, 4), mask_0f);
    let s_hi = _mm256_cvtepu8_epi16(bytes_1);
    let zmm_lo = _mm512_castsi256_si512(s_lo);
    _mm512_inserti64x4::<1>(zmm_lo, s_hi)
}

#[inline(always)]
unsafe fn buf_ring_add(z_buf_ptr: *mut __m512i, s_ptr: *const u8) {
    for i in 0..32 {
        let s_val = load_unpacked_s(s_ptr.add(i * 16));
        let z_val = _mm512_load_si512(z_buf_ptr.add(i));
        let sum = _mm512_add_epi16(z_val, s_val);
        _mm512_store_si512(z_buf_ptr.add(i), sum);
    }
}

#[inline(always)]
unsafe fn buf_ring_sub(z_buf_ptr: *mut __m512i, s_ptr: *const u8) {
    for i in 0..32 {
        let s_val = load_unpacked_s(s_ptr.add(i * 16));
        let z_val = _mm512_load_si512(z_buf_ptr.add(i));
        let diff = _mm512_sub_epi16(z_val, s_val);
        _mm512_store_si512(z_buf_ptr.add(i), diff);
    }
}

unsafe fn fold_witness_ring(z_buf: *mut __m512i, id: &[i16], s_ptr: *const u8){
    for i in 0..16 {
        let val = id[i] as usize;
        let rotate = (val >> 5) & 0x1F;
        let offset = val & 0x1F;
        let start = (offset << 6) + rotate + 1;
        let target_ptr = z_buf.add(start);
        if (val >> 10) == 0 {  // 1
            buf_ring_add(target_ptr, s_ptr);
        } else {               // -1
            buf_ring_sub(target_ptr, s_ptr);
        }
    }
}

unsafe fn fold_witness_row(z: &mut [i16], id: &[i16], s_row_ptr: *const u8, z_buf: &mut [__m512i; 2080]){
    ptr::write_bytes(z_buf.as_mut_ptr(), 0, 2080);  
    let z_buf_ptr = z_buf.as_mut_ptr();
    for i in 0..(1<<10){
        let id_slice = &id[i<<4..(i+1)<<4];
        let s_ptr_for_ring = s_row_ptr.add(i << 9); 
        fold_witness_ring(z_buf_ptr, id_slice, s_ptr_for_ring);
    }
    // Merge
    let buf_ptr = z_buf_ptr as *const i16;
    let z_ptr = z.as_mut_ptr() as *mut __m512i;
    let mut acc = [_mm512_setzero_si512(); 32]; 
    for i in 0..32 {
        let idx = (i << 11) + 32 - i;
        let ptr = buf_ptr.add(idx);
        for j in 0..32 {
            let tmp_pos = _mm512_loadu_si512(ptr.add(j << 5) as *const __m512i);
            acc[j] = _mm512_add_epi16(acc[j], tmp_pos);
            let tmp_neg = _mm512_loadu_si512(ptr.add((j << 5) + 1024) as *const __m512i);
            acc[j] = _mm512_sub_epi16(acc[j], tmp_neg);
        }
    }
    for j in 0..32 {
        _mm512_storeu_si512(z_ptr.add(j), acc[j]);
    }
}

unsafe fn fold_witness(z: &mut [i16], c: &[i16], s: &[u8]) {
    z.par_chunks_mut(1<<10).enumerate().for_each_init(|| vec![_mm512_setzero_si512(); 2080].into_boxed_slice(),|thread_local_buf, (i, z_row)| {
        let s_ptr = s.as_ptr() as *const u8;
        let s_row_ptr = s_ptr.add(i << 19);  
        let z_buf_array = &mut *(thread_local_buf.as_mut_ptr() as *mut [__m512i; 2080]);
        fold_witness_row(z_row, c, s_row_ptr, z_buf_array);
    });
}

pub unsafe fn fold_witness_u32(z: &mut [u32], c: &[i16], s: &[u8]){
    let q = ((1u64 << 32) - 99) as u32;
    let mut z_i16 = vec![0i16; z.len()];
    fold_witness(&mut z_i16, c, s);
    z.par_iter_mut().zip(z_i16.par_iter()).for_each(|(z_out, &z_in)| {
        let mut val = z_in as i64;
        if val < 0 {
            val += q as i64;
        }
        *z_out = val as u32;
    });
}