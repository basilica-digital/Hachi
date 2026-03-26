use std::ptr;

#[inline(always)]
unsafe fn buf_ring_add(z_buf_ptr: *mut __m512i, s_ptr: *const __m512i) {
    for i in 0..32 {
        let s_val = _mm512_loadu_si512(s_ptr.add(i));
        let z_val = _mm512_load_si512(z_buf_ptr.add(i));
        let sum = _mm512_add_epi16(z_val, s_val);
        _mm512_store_si512(z_buf_ptr.add(i), sum);
    }
}

#[inline(always)]
unsafe fn buf_ring_sub(z_buf_ptr: *mut __m512i, s_ptr: *const __m512i) {
    for i in 0..32 {
        let s_val = _mm512_loadu_si512(s_ptr.add(i));
        let z_val = _mm512_load_si512(z_buf_ptr.add(i));
        let diff = _mm512_sub_epi16(z_val, s_val);
        _mm512_store_si512(z_buf_ptr.add(i), diff);
    }
}

#[inline(always)]
pub unsafe fn fold_witness_ring(z_buf: *mut __m512i, id: &[i16], s_ptr: *const __m512i){
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

#[inline(always)]
pub unsafe fn fold_witness_row(z: &mut [i16], id: &[i16], s_row_ptr: *const __m512i, z_buf: &mut [__m512i; 2080]){
    ptr::write_bytes(z_buf.as_mut_ptr(), 0, 2080);  
    let z_buf_ptr = z_buf.as_mut_ptr();
    for i in 0..(1<<10){
        let id_slice = &id[i<<4..(i+1)<<4];
        let s_ptr_for_ring = s_row_ptr.add(i << 5);
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

#[inline(never)]
pub unsafe fn fold_witness(z: &mut [i16], c: &[i16], s: &[i16]) {
    z.par_chunks_mut(1<<10).enumerate().for_each_init(|| vec![_mm512_setzero_si512(); 2080].into_boxed_slice(),|thread_local_buf, (i, z_row)| {
        let s_ptr = s.as_ptr() as *const __m512i;
        let s_row_ptr = s_ptr.add(i << 15);
        let z_buf_array = &mut *(thread_local_buf.as_mut_ptr() as *mut [__m512i; 2080]);
        fold_witness_row(z_row, c, s_row_ptr, z_buf_array);
    });
}

#[inline(never)]
pub unsafe fn fold_witness_u32(z: &mut [u32], c: &[u32], s: &[u32]) {
    // convert data type
    let q = ((1u64 << 32) - 99) as u32;
    let mut c_i16 = vec![0i16; (1<<14)];
    let mut s_i16 = vec![0i16; s.len()];
    let mut z_i16 = vec![0i16; z.len()];
    for i in 0..(c.len()/1024){
        let mut id = i*16;
        for j in 0..1024{
            if c[i*1024+j] == 1{
                c_i16[id] = j as i16;
                id += 1;
            }else if c[i*1024+j] == (q - 1){
                c_i16[id] = (j+1024) as i16;
                id += 1;
            }
            if id == (i+1)*16{
                break;
            }
        }
    }
    for i in 0..s.len(){
        s_i16[i] = s[i] as i16;
    }

    // fold
    fold_witness(&mut z_i16, &c_i16, &s_i16);
    // convert data type
    for i in 0..z.len(){
        let mut val = z_i16[i] as i64;
        if val < 0 {
            val += q as i64;
        }
        z[i] = val as u32;
    }
}
