#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;
use crate::NTT::transpose::*;

pub unsafe fn intt_257(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_257(s0);
    let s2 = transpose(s1);
    let s3 = intt1234_257(s2);
    s3
}

pub unsafe fn intt1234_257(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l1 ~ l4
    // inverse Layer 4
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_257(a[2*i], a[2*i+1], l4[2*i], l4[2*i+1]);
    }
    // inverse Layer 3
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_257(a0[i*4+j], a0[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // inverse Layer 2
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_257(a1[i], a1[i+4], l2[0], l2[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_257(a1[i+8], a1[i+12], l2[2], l2[3]);
    }
    // inverse Layer 1
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_257(a0[i], a0[i+8], l1[0], l1[1]);
    }
    a1
}

pub unsafe fn intt5678_257(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // inverse Layer 8
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_257(a[2*i], a[2*i+1], l8[2*i], l8[2*i+1]);
    }
    // inverse Layer 7
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_257(a0[i*4+j], a0[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_257(a1[i], a1[i+4], l6[0], l6[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_257(a1[i+8], a1[i+12], l6[2], l6[3]);
    }
    // inverse Layer 5
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_257(a0[i], a0[i+8], l5[0], l5[1]);
    }
    a1
}

pub unsafe fn intt_3329(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_3329(s0);
    let s2 = transpose(s1);
    let s3 = intt1234_3329(s2);
    s3
}

pub unsafe fn intt1234_3329(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l1 ~ l4
    // inverse Layer 4
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_3329(a[2*i], a[2*i+1], l4[2*i], l4[2*i+1]);
    }
    // inverse Layer 3
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_3329(a0[i*4+j], a0[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // inverse Layer 2
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_3329(a1[i], a1[i+4], l2[0], l2[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_3329(a1[i+8], a1[i+12], l2[2], l2[3]);
    }
    // inverse Layer 1
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_3329(a0[i], a0[i+8], l1[0], l1[1]);
    }
    a1
}

pub unsafe fn intt5678_3329(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // inverse Layer 8
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_3329(a[2*i], a[2*i+1], l8[2*i], l8[2*i+1]);
    }
    // inverse Layer 7
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_3329(a0[i*4+j], a0[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_3329(a1[i], a1[i+4], l6[0], l6[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_3329(a1[i+8], a1[i+12], l6[2], l6[3]);
    }
    // inverse Layer 5
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_3329(a0[i], a0[i+8], l5[0], l5[1]);
    }
    a1
}

pub unsafe fn intt_7681(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_7681(s0);
    let s2 = transpose(s1);
    let s3 = intt1234_7681(s2);
    s3
}

pub unsafe fn intt1234_7681(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l1 ~ l4
    // inverse Layer 4
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_7681(a[2*i], a[2*i+1], l4[2*i], l4[2*i+1]);
    }
    // inverse Layer 3
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_7681(a0[i*4+j], a0[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // inverse Layer 2
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_7681(a1[i], a1[i+4], l2[0], l2[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_7681(a1[i+8], a1[i+12], l2[2], l2[3]);
    }
    // inverse Layer 1
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_7681(a0[i], a0[i+8], l1[0], l1[1]);
    }
    a1
}

pub unsafe fn intt5678_7681(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // inverse Layer 8
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_7681(a[2*i], a[2*i+1], l8[2*i], l8[2*i+1]);
    }
    // inverse Layer 7
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_7681(a0[i*4+j], a0[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_7681(a1[i], a1[i+4], l6[0], l6[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_7681(a1[i+8], a1[i+12], l6[2], l6[3]);
    }
    // inverse Layer 5
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_7681(a0[i], a0[i+8], l5[0], l5[1]);
    }
    a1
}

pub unsafe fn intt_7937(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_7937(s0);
    let s2 = transpose(s1);
    let s3 = intt1234_7937(s2);
    s3
}

pub unsafe fn intt1234_7937(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l1 ~ l4
    // inverse Layer 4
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_7937(a[2*i], a[2*i+1], l4[2*i], l4[2*i+1]);
    }
    // inverse Layer 3
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_7937(a0[i*4+j], a0[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // inverse Layer 2
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_7937(a1[i], a1[i+4], l2[0], l2[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_7937(a1[i+8], a1[i+12], l2[2], l2[3]);
    }
    // inverse Layer 1
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_7937(a0[i], a0[i+8], l1[0], l1[1]);
    }
    a1
}

pub unsafe fn intt5678_7937(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // inverse Layer 8
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_7937(a[2*i], a[2*i+1], l8[2*i], l8[2*i+1]);
    }
    // inverse Layer 7
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_7937(a0[i*4+j], a0[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_7937(a1[i], a1[i+4], l6[0], l6[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_7937(a1[i+8], a1[i+12], l6[2], l6[3]);
    }
    // inverse Layer 5
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_7937(a0[i], a0[i+8], l5[0], l5[1]);
    }
    a1
}

pub unsafe fn intt_9473(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_9473(s0);
    let s2 = transpose(s1);
    let s3 = intt1234_9473(s2);
    s3
}

pub unsafe fn intt1234_9473(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l1 ~ l4
    // inverse Layer 4
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_9473(a[2*i], a[2*i+1], l4[2*i], l4[2*i+1]);
    }
    // inverse Layer 3
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_9473(a0[i*4+j], a0[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // inverse Layer 2
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_9473(a1[i], a1[i+4], l2[0], l2[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_9473(a1[i+8], a1[i+12], l2[2], l2[3]);
    }
    // inverse Layer 1
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_9473(a0[i], a0[i+8], l1[0], l1[1]);
    }
    a1
}

pub unsafe fn intt5678_9473(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // inverse Layer 8
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_9473(a[2*i], a[2*i+1], l8[2*i], l8[2*i+1]);
    }
    // inverse Layer 7
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_9473(a0[i*4+j], a0[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_9473(a1[i], a1[i+4], l6[0], l6[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_9473(a1[i+8], a1[i+12], l6[2], l6[3]);
    }
    // inverse Layer 5
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_9473(a0[i], a0[i+8], l5[0], l5[1]);
    }
    a1
}

pub unsafe fn intt_10753(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_10753(s0);
    let s2 = transpose(s1);
    let s3 = intt1234_10753(s2);
    s3
}

pub unsafe fn intt1234_10753(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l1 ~ l4
    // inverse Layer 4
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_10753(a[2*i], a[2*i+1], l4[2*i], l4[2*i+1]);
    }
    // inverse Layer 3
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_10753(a0[i*4+j], a0[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // inverse Layer 2
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_10753(a1[i], a1[i+4], l2[0], l2[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_10753(a1[i+8], a1[i+12], l2[2], l2[3]);
    }
    // inverse Layer 1
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_10753(a0[i], a0[i+8], l1[0], l1[1]);
    }
    a1
}

pub unsafe fn intt5678_10753(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // inverse Layer 8
    for i in 0..8{
        [a0[2*i], a0[2*i+1]] = barrett_ibutterfly_10753(a[2*i], a[2*i+1], l8[2*i], l8[2*i+1]);
    }
    // inverse Layer 7
    for i in 0..4{
        for j in 0..2{
            [a1[i*4+j], a1[i*4+j+2]] = barrett_butterfly_10753(a0[i*4+j], a0[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        [a0[i], a0[i+4]] = barrett_butterfly_10753(a1[i], a1[i+4], l6[0], l6[1]);
        [a0[i+8], a0[i+12]] = barrett_butterfly_10753(a1[i+8], a1[i+12], l6[2], l6[3]);
    }
    // inverse Layer 5
    for i in 0..8{
        [a1[i], a1[i+8]] = barrett_butterfly_10753(a0[i], a0[i+8], l5[0], l5[1]);
    }
    a1
}