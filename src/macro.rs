use std::ops::{Deref, DerefMut};

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

pub struct AlignedU32Vec {
    inner: Vec<Align64<[u32; 16]>>,
    len: usize, 
}

impl Deref for AlignedU32Vec {
    type Target = [u32];

    fn deref(&self) -> &Self::Target {
        unsafe {
            std::slice::from_raw_parts(self.inner.as_ptr() as *const u32, self.len)
        }
    }
}

impl DerefMut for AlignedU32Vec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe {
            std::slice::from_raw_parts_mut(self.inner.as_mut_ptr() as *mut u32, self.len)
        }
    }
}

pub fn generate_random_data_4bit(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    
    data.par_iter_mut().for_each(|chunk| {
        let mut rng = rand::thread_rng();
        for x in chunk.0.iter_mut() {
            *x = rng.next_u32() & 0xF;
        }
    });
    
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn generate_random_data_2bit(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    
    data.par_iter_mut().for_each(|chunk| {
        let mut rng = rand::thread_rng();
        for x in chunk.0.iter_mut() {
            *x = rng.next_u32() & 0x3;
        }
    });
    
    AlignedU32Vec { inner: data, len: total_u32 }
}


pub fn generate_random_data_32bit(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    
    data.par_iter_mut().for_each(|chunk| {
        let mut rng = rand::thread_rng();
        for x in chunk.0.iter_mut() {
            *x = rng.next_u32() & 0xFFFFFFFF;
        }
    });
    
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn generate_random_data_30bit(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    
    data.par_iter_mut().for_each(|chunk| {
        let mut rng = rand::thread_rng();
        for x in chunk.0.iter_mut() {
            *x = rng.next_u32() & 0x3FFFFFFF;
        }
    });
    
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn generate_random_eval_points_30bit(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    let ring_size = n / 16;
    
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    
    //data.par_iter_mut().for_each(|chunk| {
    for i in 0..sz {
        let mut rng = rand::thread_rng();
        data[i*ring_size].0[0] = rng.next_u32() & 0x3FFFFFFF;
    }
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn sparse_random_1(s: &mut [u64]){
    let mut rng = thread_rng();
    // 2 + 16 + 10
    for i in 0..(1<<28){
        let index: usize = rng.gen_range(0..256);
        s[i*4+index/64] = 1u64<<(index%64);
    }
}

pub fn generate_sparse_c(num_rings: usize) -> Vec<u32> {
    let n = 1024;
    let weight = 16;
    let q = (1u64 << 32) - 99;
    
    let mut c = vec![0u32; num_rings * n];
    let mut rng = rand::thread_rng();

    for i in 0..num_rings {
        let indices = rand::seq::index::sample(&mut rng, n, weight);
        for idx in indices.into_iter() {
            let sign = rng.gen_bool(0.5);
            c[i*n+idx] = if sign { 1 } else { (q - 1) as u32 };
        }
    }
    c
}

pub fn generate_sparse_c_idx(num_rings: usize) -> Vec<i16> {
    let n = 16;
    let mut c = vec![0i16; num_rings*n];
    let mut rng = rand::thread_rng();

    for i in 0..num_rings*n{
        c[i] = rng.gen_range(0..2048 as i16);
    }
    c
}


pub fn sparse_random_1_index(s: &mut [u8]){
    let mut rng = thread_rng();
    for i in 0..s.len(){
        let index: usize = rng.gen_range(0..256);
        s[i] = index as u8;
    }
}