fn mod_inverse(a: u64, m: u64) -> u64 {
    let mut mn = (m as i64, a as i64);
    let mut xy = (0, 1);
    while mn.1 != 0 {
        xy = (xy.1, xy.0 - (mn.0 / mn.1) * xy.1);
        mn = (mn.1, mn.0 % mn.1);
    }
    if xy.0 < 0 {
        (xy.0 + m as i64) as u64
    } else {
        xy.0 as u64
    }
}

#[inline(always)]
pub unsafe fn load_vec_to_m512(src: &[u32], offset: usize) -> __m512i {
    let ptr = src.as_ptr().add(offset) as *const __m512i;
    _mm512_load_si512(ptr)
}

#[inline(always)]
pub unsafe fn store_m512_to_vec(src: __m512i, dst: &mut [u32], offset: usize) {
    let ptr = dst.as_mut_ptr().add(offset) as *mut __m512i;
    _mm512_store_si512(ptr, src)
}

pub fn generate_random_data_4bit(sz: usize, n: usize) -> Vec<u32> {
    let total = sz * n;
    let mut data = vec![0u32; total];
    data.par_iter_mut().for_each(|x| {
        let mut rng = thread_rng();
        *x = rng.next_u32() & 0xF;
    });
    data
}

pub fn generate_random_data_30bit(sz: usize, n: usize) -> Vec<u32> {
    let total = sz * n;
    let mut data = vec![0u32; total];
    data.par_iter_mut().for_each(|x| {
        let mut rng = thread_rng();
        *x = rng.next_u32() & 0x3FFFFFFF; 
    });
    data
}

pub fn sparse_random_1(s: &mut [u64]){
    let mut rng = thread_rng();
    // 2 + 16 + 10
    for i in 0..(1<<28){
        let index: usize = rng.gen_range(0..256);
        s[i*4+index/64] = 1u64<<(index%64);
    }
}


pub fn sparse_random_1_index(s: &mut [u8]){
    let mut rng = thread_rng();
    for i in 0..s.len(){
        let index: usize = rng.gen_range(0..256);
        s[i] = index as u8;
    }
}