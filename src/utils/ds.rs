use std::ops::{Deref, DerefMut};
use rand::seq::SliceRandom;

pub struct AlignedI8Vec {
    pub inner: Vec<Align64<[i8; 64]>>,
    pub len: usize, 
}

impl std::ops::Deref for AlignedI8Vec {
    type Target = [i8];
    fn deref(&self) -> &Self::Target {
        unsafe { std::slice::from_raw_parts(self.inner.as_ptr() as *const i8, self.len) }
    }
}

impl std::ops::DerefMut for AlignedI8Vec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe { std::slice::from_raw_parts_mut(self.inner.as_mut_ptr() as *mut i8, self.len) }
    }
}

pub struct AlignedU8Vec {
    pub inner: Vec<Align64<[u8; 64]>>,
    pub len: usize, 
}

impl std::ops::Deref for AlignedU8Vec {
    type Target = [u8];
    fn deref(&self) -> &Self::Target {
        unsafe { std::slice::from_raw_parts(self.inner.as_ptr() as *const u8, self.len) }
    }
}

impl std::ops::DerefMut for AlignedU8Vec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe { std::slice::from_raw_parts_mut(self.inner.as_mut_ptr() as *mut u8, self.len) }
    }
}


#[repr(C, align(64))]
#[derive(Clone, Copy)]
pub struct Align64<T>(pub T);

#[repr(align(64))]
#[derive(Clone)]
pub struct Align64U32(pub [u32; 16]); // 16 * 4 bytes = 64 bytes

pub struct AlignedU32Vec {
    pub inner: Vec<Align64<[u32; 16]>>,
    pub len: usize, 
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