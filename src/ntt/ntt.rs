#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;
use crate::NTT::transpose::*;

pub unsafe fn ntt_257(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_257(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_257(s2);
    s3
}

pub unsafe fn ntt1234_257(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let l0:[__m256i;2] = [_mm256_set1_epi16(16), _mm256_set1_epi16(4080)];
    let l1:[__m256i;4] = [_mm256_set1_epi16(4), _mm256_set1_epi16(1020), _mm256_set1_epi16(64), _mm256_set1_epi16(16320)];
    let l2:[__m256i;8] = [_mm256_set1_epi16(2), _mm256_set1_epi16(510), _mm256_set1_epi16(32), _mm256_set1_epi16(8160), _mm256_set1_epi16(8), _mm256_set1_epi16(2040), _mm256_set1_epi16(128), _mm256_set1_epi16(32640)];
    let l3:[__m256i;16] = [_mm256_set1_epi16(60), _mm256_set1_epi16(15300), _mm256_set1_epi16(68), _mm256_set1_epi16(17340), _mm256_set1_epi16(17), _mm256_set1_epi16(4335), _mm256_set1_epi16(15), _mm256_set1_epi16(3825), _mm256_set1_epi16(120), _mm256_set1_epi16(30600), _mm256_set1_epi16(121), _mm256_set1_epi16(30855), _mm256_set1_epi16(34), _mm256_set1_epi16(8670), _mm256_set1_epi16(30), _mm256_set1_epi16(7650)];

    // Layer 1
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_257(a[i], a[i+8], l1[0], l1[1]);
    }
    // Layer 2
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_257(a0[i], a0[i+4], l2[0], l2[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_257(a0[i+8], a0[i+12], l2[2], l2[3]);
    }
    // Layer 3
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_257(a1[i*4+j], a1[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // Layer 4
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_257(a0[i*2], a0[i*2+1], l4[2*i], l4[2*i+1]);
    }
    a1
}

pub unsafe fn ntt5678_257(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // Layer 5
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_257(a[i], a[i+8], l5[0], l5[1]);
    }
    // Layer 6
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_257(a0[i], a0[i+4], l6[0], l6[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_257(a0[i+8], a0[i+12], l6[2], l6[3]);
    }
    // Layer 7
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_257(a1[i*4+j], a1[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // Layer 8
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_257(a0[i*2], a0[i*2+1], l8[2*i], l8[2*i+1]);
    }
    a1
}

pub unsafe fn ntt_3329(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_3329(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_3329(s2);
    s3
}

pub unsafe fn ntt1234_3329(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let l0:[__m256i;2] = [_mm256_set1_epi16(1600), _mm256_set1_epi16(31498)];
    let l1:[__m256i;4] = [_mm256_set1_epi16(40), _mm256_set1_epi16(787), _mm256_set1_epi16(749), _mm256_set1_epi16(14745)];
    let l2:[__m256i;8] = [_mm256_set1_epi16(848), _mm256_set1_epi16(16694), _mm256_set1_epi16(1432), _mm256_set1_epi16(28191), _mm256_set1_epi16(630), _mm256_set1_epi16(12402), _mm256_set1_epi16(687), _mm256_set1_epi16(13525)];
    let l3:[__m256i;16] = [_mm256_set1_epi16(569), _mm256_set1_epi16(11202), _mm256_set1_epi16(1583), _mm256_set1_epi16(31164), _mm256_set1_epi16(69), _mm256_set1_epi16(1358), _mm256_set1_epi16(543), _mm256_set1_epi16(10690), _mm256_set1_epi16(193), _mm256_set1_epi16(3799), _mm256_set1_epi16(797), _mm256_set1_epi16(15690), _mm256_set1_epi16(1410), _mm256_set1_epi16(27758), _mm256_set1_epi16(1062), _mm256_set1_epi16(20907)];

    // Layer 1
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_3329(a[i], a[i+8], l1[0], l1[1]);
    }
    // Layer 2
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_3329(a0[i], a0[i+4], l2[0], l2[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_3329(a0[i+8], a0[i+12], l2[2], l2[3]);
    }
    // Layer 3
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_3329(a1[i*4+j], a1[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // Layer 4
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_3329(a0[i*2], a0[i*2+1], l4[2*i], l4[2*i+1]);
    }
    a1
}

pub unsafe fn ntt5678_3329(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // Layer 5
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_3329(a[i], a[i+8], l5[0], l5[1]);
    }
    // Layer 6
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_3329(a0[i], a0[i+4], l6[0], l6[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_3329(a0[i+8], a0[i+12], l6[2], l6[3]);
    }
    // Layer 7
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_3329(a1[i*4+j], a1[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // Layer 8
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_3329(a0[i*2], a0[i*2+1], l8[2*i], l8[2*i+1]);
    }
    a1
}

pub unsafe fn ntt_7681(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_7681(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_7681(s2);
    s3
}

pub unsafe fn ntt1234_7681(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let l0:[__m256i;2] = [_mm256_set1_epi16(3383), _mm256_set1_epi16(28865)];
    let l1:[__m256i;4] = [_mm256_set1_epi16(1925), _mm256_set1_epi16(16425), _mm256_set1_epi16(1213), _mm256_set1_epi16(10350)];
    let l2:[__m256i;8] = [_mm256_set1_epi16(583), _mm256_set1_epi16(4974), _mm256_set1_epi16(1728), _mm256_set1_epi16(14744), _mm256_set1_epi16(527), _mm256_set1_epi16(4496), _mm256_set1_epi16(849), _mm256_set1_epi16(7244)];
    let l3:[__m256i;16] = [_mm256_set1_epi16(2381), _mm256_set1_epi16(20315), _mm256_set1_epi16(2446), _mm256_set1_epi16(20870), _mm256_set1_epi16(97), _mm256_set1_epi16(828), _mm256_set1_epi16(2132), _mm256_set1_epi16(18191), _mm256_set1_epi16(2784), _mm256_set1_epi16(23754), _mm256_set1_epi16(1366), _mm256_set1_epi16(11655), _mm256_set1_epi16(2138), _mm256_set1_epi16(18242), _mm256_set1_epi16(2648), _mm256_set1_epi16(22593)];

    // Layer 1
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_7681(a[i], a[i+8], l1[0], l1[1]);
    }
    // Layer 2
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_7681(a0[i], a0[i+4], l2[0], l2[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_7681(a0[i+8], a0[i+12], l2[2], l2[3]);
    }
    // Layer 3
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_7681(a1[i*4+j], a1[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // Layer 4
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_7681(a0[i*2], a0[i*2+1], l4[2*i], l4[2*i+1]);
    }
    a1
}

pub unsafe fn ntt5678_7681(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // Layer 5
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_7681(a[i], a[i+8], l5[0], l5[1]);
    }
    // Layer 6
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_7681(a0[i], a0[i+4], l6[0], l6[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_7681(a0[i+8], a0[i+12], l6[2], l6[3]);
    }
    // Layer 7
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_7681(a1[i*4+j], a1[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // Layer 8
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_7681(a0[i*2], a0[i*2+1], l8[2*i], l8[2*i+1]);
    }
    a1
}

pub unsafe fn ntt_7937(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_7937(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_7937(s2);
    s3
}

pub unsafe fn ntt1234_7937(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let l0:[__m256i;2] = [_mm256_set1_epi16(1962), _mm256_set1_epi16(16200)];
    let l1:[__m256i;4] = [_mm256_set1_epi16(3323), _mm256_set1_epi16(27438), _mm256_set1_epi16(3449), _mm256_set1_epi16(28478)];
    let l2:[__m256i;8] = [_mm256_set1_epi16(2781), _mm256_set1_epi16(22963), _mm256_set1_epi16(3603), _mm256_set1_epi16(29750), _mm256_set1_epi16(2595), _mm256_set1_epi16(21427), _mm256_set1_epi16(3773), _mm256_set1_epi16(31154)];
    let l3:[__m256i;16] = [_mm256_set1_epi16(888), _mm256_set1_epi16(7332), _mm256_set1_epi16(3884), _mm256_set1_epi16(32070), _mm256_set1_epi16(1740), _mm256_set1_epi16(14367), _mm256_set1_epi16(970), _mm256_set1_epi16(8009), _mm256_set1_epi16(1121), _mm256_set1_epi16(9256), _mm256_set1_epi16(853), _mm256_set1_epi16(7043), _mm256_set1_epi16(2630), _mm256_set1_epi16(21716), _mm256_set1_epi16(1010), _mm256_set1_epi16(8340)];

    // Layer 1
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_7937(a[i], a[i+8], l1[0], l1[1]);
    }
    // Layer 2
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_7937(a0[i], a0[i+4], l2[0], l2[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_7937(a0[i+8], a0[i+12], l2[2], l2[3]);
    }
    // Layer 3
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_7937(a1[i*4+j], a1[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // Layer 4
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_7937(a0[i*2], a0[i*2+1], l4[2*i], l4[2*i+1]);
    }
    a1
}

pub unsafe fn ntt5678_7937(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // Layer 5
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_7937(a[i], a[i+8], l5[0], l5[1]);
    }
    // Layer 6
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_7937(a0[i], a0[i+4], l6[0], l6[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_7937(a0[i+8], a0[i+12], l6[2], l6[3]);
    }
    // Layer 7
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_7937(a1[i*4+j], a1[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // Layer 8
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_7937(a0[i*2], a0[i*2+1], l8[2*i], l8[2*i+1]);
    }
    a1
}

pub unsafe fn ntt_9473(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_9473(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_9473(s2);
    s3
}

pub unsafe fn ntt1234_9473(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let l0:[__m256i;2] = [_mm256_set1_epi16(1172), _mm256_set1_epi16(8108)];
    let l1:[__m256i;4] = [_mm256_set1_epi16(4677), _mm256_set1_epi16(32356), _mm256_set1_epi16(3423), _mm256_set1_epi16(23681)];
    let l2:[__m256i;8] = [_mm256_set1_epi16(1906), _mm256_set1_epi16(13186), _mm256_set1_epi16(1796), _mm256_set1_epi16(12425), _mm256_set1_epi16(2659), _mm256_set1_epi16(18395), _mm256_set1_epi16(269), _mm256_set1_epi16(1861)];
    let l3:[__m256i;16] = [_mm256_set1_epi16(3221), _mm256_set1_epi16(22283), _mm256_set1_epi16(4715), _mm256_set1_epi16(32619), _mm256_set1_epi16(1089), _mm256_set1_epi16(7534), _mm256_set1_epi16(2547), _mm256_set1_epi16(17621), _mm256_set1_epi16(4406), _mm256_set1_epi16(30482), _mm256_set1_epi16(1047), _mm256_set1_epi16(7243), _mm256_set1_epi16(722), _mm256_set1_epi16(4995), _mm256_set1_epi16(3087), _mm256_set1_epi16(21356)];


    // Layer 1
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_9473(a[i], a[i+8], l1[0], l1[1]);
    }
    // Layer 2
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_9473(a0[i], a0[i+4], l2[0], l2[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_9473(a0[i+8], a0[i+12], l2[2], l2[3]);
    }
    // Layer 3
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_9473(a1[i*4+j], a1[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // Layer 4
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_9473(a0[i*2], a0[i*2+1], l4[2*i], l4[2*i+1]);
    }
    a1
}

pub unsafe fn ntt5678_9473(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // Layer 5
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_9473(a[i], a[i+8], l5[0], l5[1]);
    }
    // Layer 6
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_9473(a0[i], a0[i+4], l6[0], l6[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_9473(a0[i+8], a0[i+12], l6[2], l6[3]);
    }
    // Layer 7
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_9473(a1[i*4+j], a1[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // Layer 8
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_9473(a0[i*2], a0[i*2+1], l8[2*i], l8[2*i+1]);
    }
    a1
}

pub unsafe fn ntt_10753(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_10753(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_10753(s2);
    s3
}

pub unsafe fn ntt1234_10753(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let l0:[__m256i;2] = [_mm256_set1_epi16(4489), _mm256_set1_epi16(27359)];
    let l1:[__m256i;4] = [_mm256_set1_epi16(67), _mm256_set1_epi16(408), _mm256_set1_epi16(321), _mm256_set1_epi16(1956)];
    let l2:[__m256i;8] = [_mm256_set1_epi16(3422), _mm256_set1_epi16(20856), _mm256_set1_epi16(4679), _mm256_set1_epi16(28517), _mm256_set1_epi16(1656), _mm256_set1_epi16(10093), _mm256_set1_epi16(3461), _mm256_set1_epi16(21094)];
    let l3:[__m256i;16] = [_mm256_set1_epi16(1560), _mm256_set1_epi16(9508), _mm256_set1_epi16(2637), _mm256_set1_epi16(16072), _mm256_set1_epi16(4631), _mm256_set1_epi16(28224), _mm256_set1_epi16(3010), _mm256_set1_epi16(18345), _mm256_set1_epi16(2640), _mm256_set1_epi16(16090), _mm256_set1_epi16(1154), _mm256_set1_epi16(7033), _mm256_set1_epi16(4832), _mm256_set1_epi16(29449), _mm256_set1_epi16(2047), _mm256_set1_epi16(12476)];


    // Layer 1
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_10753(a[i], a[i+8], l1[0], l1[1]);
    }
    // Layer 2
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_10753(a0[i], a0[i+4], l2[0], l2[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_10753(a0[i+8], a0[i+12], l2[2], l2[3]);
    }
    // Layer 3
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_10753(a1[i*4+j], a1[i*4+j+2], l3[2*i], l3[2*i+1]);
        }
    }
    // Layer 4
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_10753(a0[i*2], a0[i*2+1], l4[2*i], l4[2*i+1]);
    }
    a1
}

pub unsafe fn ntt5678_10753(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    // let l5 ~ l8
    // Layer 5
    for i in 0..8{
        [a0[i], a0[i+8]] = barrett_butterfly_10753(a[i], a[i+8], l5[0], l5[1]);
    }
    // Layer 6
    for i in 0..4{
        [a1[i], a1[i+4]] = barrett_butterfly_10753(a0[i], a0[i+4], l6[0], l6[1]);
        [a1[i+8], a1[i+12]] = barrett_butterfly_10753(a0[i+8], a0[i+12], l6[2], l6[3]);
    }
    // Layer 7
    for i in 0..4{
        for j in 0..2{
            [a0[i*4+j], a0[i*4+j+2]] = barrett_butterfly_10753(a1[i*4+j], a1[i*4+j+2], l7[2*i], l7[2*i+1]);
        }
    }
    // Layer 8
    for i in 0..8{
        [a1[i*2], a1[i*2+1]] = barrett_butterfly_10753(a0[i*2], a0[i*2+1], l8[2*i], l8[2*i+1]);
    }
    a1
}