#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

mod field;
use crate::field::fields::*;
use crate::field::NTT::transpose;
// use std::mem;

pub unsafe fn ntt_26113(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = ntt1234_26113(s0);
    let s2 = transpose(s1);
    let s3 = ntt5678_26113(s2);
    s3
}

pub unsafe fn intt_26113(s0: [__m256i;16]) -> [__m256i;16]{
    let s1 = intt5678_26113(s0);
    let s2 = transpose(s1);
    let s3 = intt1324_26113(s2);
    s3
}



// twiddle factor: [8534,, 7720, 619,, 5428, 1910, 8645, 7205,, 5700, 4719, 3045, 3595, 7249, 1269, 2121, 4305] 
// 9604, 8171, 8888, 8073, 12760, 2630, 8964, 12314, 8341, 1944, 2138, 7295, 4994, 2380, 10892, 9952, 98, 716, 8436, 717, 1039, 11594, 9684, 4389, 6203, 5351, 4082, 1046, 10227, 7572, 11167, 12841, 3217, 9115, 6735, 1777, 7915, 7721, 9859, 620, 9370, 5574, 3390, 2944, 7764, 9295, 1124, 8745, 6147, 2519, 7519, 7505, 860, 1487, 10080, 6498, 7400, 10366, 10825, 7244, 5746, 3850, 5406, 6867, 6344, 7447, 12308, 9986, 608, 7815, 6580, 10770, 5705, 11838, 10031, 6140, 7732, 2663, 3278, 7429, 1121, 9256, 10717, 11152, 3122, 7888, 459, 156, 7358, 8593, 10940, 7985, 5000, 1358, 4986, 12447, 11718, 11378, 7528, 5972, 5964, 2539, 9783, 4861, 10092, 4454, 10952, 5941, 1807, 11845, 5698, 4326, 6077, 700, 1391, 10621, 5237, 12898, 6716, 3691, 458, 8378, 10505, 3741, 9766, 9652, 5289, 13051, 5232, 3342, 5771, 596, 2924, 10612, 8159, 11648, 11759, 953, 6706, 10692, 1354, 13023, 2510, 7680, 9921, 7468, 691, 4544, 8928, 6182, 11953, 9524, 3296, 4363, 3410, 11058, 2127, 3283, 4637, 10963, 9739, 5053, 5753, 3662, 9034, 10580, 5343, 3864, 7342, 11541, 1036, 11083, 539, 3938, 9113, 5828, 6733, 10822, 12223, 10353, 11476, 12434, 6689, 908, 2320, 5226, 3118, 135, 8010, 6494, 1616, 3280, 3176, 1310, 7469, 1387, 11757, 8092, 4748, 7944, 858, 10532, 8942, 8842, 9110, 6339, 6891, 1318, 10645, 2697, 8779, 1789, 7009, 10077, 3344, 3813, 7657, 10112, 7792, 12883, 1690, 8084, 9700, 1590, 798, 5361, 2108, 2185, 9626, 3214, 4878, 4730, 12384, 5745, 11526, 4787, 5490, 4938, 1401, 3620]



pub unsafe fn ntt1234_26113(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut tmp: __m256i = _mm256_setzero_si256();
    // Layer 1
    for i in 0..8{
        tmp = barrett_mul_8534_26113(a[i+8]);
        a0[i] = barrett_fake_26113(_mm256_add_epi16(a[i], tmp));
        a0[i+8] = barrett_fake_26113(_mm256_sub_epi16(a[i], tmp));
    }
    // Layer 2
    for i in 0..4{
        tmp = barrett_mul_7720_26113(a0[i+4]);
        a1[i] = barrett_fake_26113(_mm256_add_epi16(a0[i], tmp));
        a1[i+4] = barrett_fake_26113(_mm256_sub_epi16(a0[i], tmp));
        tmp = barrett_mul_619_26113(a0[i+12]);
        a1[i+8] = barrett_fake_26113(_mm256_add_epi16(a0[i+8], tmp));
        a1[i+12] = barrett_fake_26113(_mm256_sub_epi16(a0[i+8], tmp));
    }
    // Layer 3
    for i in 0..2{
        tmp = barrett_mul_5428_26113(a1[i+2]);
        a0[i] = barrett_fake_26113(_mm256_add_epi16(a1[i], tmp));
        a0[i+2] = barrett_fake_26113(_mm256_sub_epi16(a1[i], tmp));
        tmp = barrett_mul_1910_26113(a1[i+6]);
        a0[i+4] = barrett_fake_26113(_mm256_add_epi16(a1[+4], tmp));
        a0[i+6] = barrett_fake_26113(_mm256_sub_epi16(a1[i+4], tmp));
        tmp = barrett_mul_8645_26113(a1[i+10]);
        a0[i+8] = barrett_fake_26113(_mm256_add_epi16(a1[i+8], tmp));
        a0[i+10] = barrett_fake_26113(_mm256_sub_epi16(a1[i+8], tmp));
        tmp = barrett_mul_7205_26113(a1[i+14]);
        a0[i+12] = barrett_fake_26113(_mm256_add_epi16(a1[i+12], tmp));
        a0[i+14] = barrett_fake_26113(_mm256_sub_epi16(a1[i+12], tmp));
    }
    // Layer 4
    tmp = barrett_mul_5700_26113(a0[1]);
    a1[0] = barrett_fake_26113(_mm256_add_epi16(a0[0], tmp));
    a1[1] = barrett_fake_26113(_mm256_sub_epi16(a0[0], tmp));
    tmp = barrett_mul_4719_26113(a0[3]);
    a1[2] = barrett_fake_26113(_mm256_add_epi16(a0[2], tmp));
    a1[3] = barrett_fake_26113(_mm256_sub_epi16(a0[2], tmp));
    tmp = barrett_mul_3045_26113(a0[5]);
    a1[4] = barrett_fake_26113(_mm256_add_epi16(a0[4], tmp));
    a1[5] = barrett_fake_26113(_mm256_sub_epi16(a0[4], tmp));
    tmp = barrett_mul_3595_26113(a0[7]);
    a1[6] = barrett_fake_26113(_mm256_add_epi16(a0[6], tmp));
    a1[7] = barrett_fake_26113(_mm256_sub_epi16(a0[6], tmp));
    tmp = barrett_mul_7249_26113(a0[9]);
    a1[8] = barrett_fake_26113(_mm256_add_epi16(a0[8], tmp));
    a1[9] = barrett_fake_26113(_mm256_sub_epi16(a0[8], tmp));
    tmp = barrett_mul_1269_26113(a0[11]);
    a1[10] = barrett_fake_26113(_mm256_add_epi16(a0[10], tmp));
    a1[11] = barrett_fake_26113(_mm256_sub_epi16(a0[10], tmp));
    tmp = barrett_mul_2121_26113(a0[13]);
    a1[12] = barrett_fake_26113(_mm256_add_epi16(a0[12], tmp));
    a1[13] = barrett_fake_26113(_mm256_sub_epi16(a0[12], tmp));
    tmp = barrett_mul_4305_26113(a0[15]);
    a1[14] = barrett_fake_26113(_mm256_add_epi16(a0[14], tmp));
    a1[15] = barrett_fake_26113(_mm256_sub_epi16(a0[14], tmp));
    a1
}

pub unsafe fn ntt5678_26113(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut tmp: __m256i = _mm256_setzero_si256();
    let twiddle5: __m256i = _mm256_set_epi16(9952, 10892, 2380, 4994, 7295, 2138, 1944, 8341, 12314, 8964, 2630, 12760, 8073, 8888, 8171, 9604);
    let twiddle5_r_overq: __m256i = _mm256_set_epi16(24976, 27335, 5973, 12533, 18308, 5365, 4878, 20933, 30904, 22497, 6600, 32023, 20260, 22306, 20506, 24103);
    let twiddle6: [__m256i;2] = [_mm256_set_epi16(1124, 7764, 3390, 9370, 9859, 7915, 6735, 3217, 11167, 10227, 4082, 6203, 9684, 1039, 8436, 98), _mm256_set_epi16(8745, 9295, 2944, 5574, 620, 7721, 1777, 9115, 12841, 7572, 1046, 5351, 4389, 11594, 717, 716)];
    let twiddle6_r_overq: [__m256i;2] = [_mm256_set_epi16(2820, 19485, 8507, 23515, 24743, 19864, 16902, 8073, 28025, 25666, 10244, 15567, 24304, 2607, 21171, 245), _mm256_set_epi16(21947, 23327, 7388, 13989, 1556, 19377, 4459, 22875, 32227, 19003, 2625, 13429, 11015, 29097, 1799, 1796)];
    let twiddle7: [__m256i;4] = [_mm256_set_epi16(6867, 5406, 3850, 5746, 7244, 10825, 10366, 7400, 6498, 10080, 1487, 860, 7505, 7519, 2519, 6147), _mm256_set_epi16(9986, 12308, 7447, 6344, 6867, 5406, 3850, 5746, 7244, 10825, 10366, 7400, 6498, 10080, 1487, 860), _mm256_set_epi16(10770, 6580, 7815, 608, 9986, 12308, 7447, 6344, 6867, 5406, 3850, 5746, 7244, 10825, 10366, 7400), _mm256_set_epi16(6140, 10031, 11838, 5705, 10770, 6580, 7815, 608, 9986, 12308, 7447, 6344, 6867, 5406, 3850, 5746)];
    let twiddle7_r_overq: [__m256i;4] = [_mm256_set_epi16(17234, 13567, 9662, 14420, 18180, 27167, 26015, 18571, 16308, 25297, 3731, 2158, 18835, 18870, 6321, 15427), _mm256_set_epi16(25061, 30889, 18689, 15921, 17234, 13567, 9662, 14420, 18180, 27167, 26015, 18571, 16308, 25297, 3731, 2158), _mm256_set_epi16(27029, 16513, 19613, 1525, 25061, 30889, 18689, 15921, 17234, 13567, 9662, 14420, 18180, 27167, 26015, 18571), _mm256_set_epi16(15409, 25174, 29709, 14317, 27029, 16513, 19613, 1525, 25061, 30889, 18689, 15921, 17234, 13567, 9662, 14420)];
    let twiddle8: [__m256i;8] = [_mm256_set_epi16(13051, 5289, 9652, 9766, 3741, 10505, 8378, 458, 3691, 6716, 12898, 5237, 10621, 1391, 700, 6077), _mm256_set_epi16(11648, 8159, 10612, 2924, 596, 5771, 3342, 5232, 13051, 5289, 9652, 9766, 3741, 10505, 8378, 458), _mm256_set_epi16(7680, 2510, 13023, 1354, 10692, 6706, 953, 11759, 11648, 8159, 10612, 2924, 596, 5771, 3342, 5232), _mm256_set_epi16(9524, 11953, 6182, 8928, 4544, 691, 7468, 9921, 7680, 2510, 13023, 1354, 10692, 6706, 953, 11759), _mm256_set_epi16(10963, 4637, 3283, 2127, 11058, 3410, 4363, 3296, 9524, 11953, 6182, 8928, 4544, 691, 7468, 9921), _mm256_set_epi16(3864, 5343, 10580, 9034, 3662, 5753, 5053, 9739, 10963, 4637, 3283, 2127, 11058, 3410, 4363, 3296), _mm256_set_epi16(5828, 9113, 3938, 539, 11083, 1036, 11541, 7342, 3864, 5343, 10580, 9034, 3662, 5753, 5053, 9739), _mm256_set_epi16(908, 6689, 12434, 11476, 10353, 12223, 10822, 6733, 5828, 9113, 3938, 539, 11083, 1036, 11541, 7342)];
    let twiddle8_r_overq: [__m256i;8] = [_mm256_set_epi16(32754, 13273, 24223, 24509, 9388, 26364, 21026, 1149, 9263, 16855, 32370, 13143, 26655, 3491, 1756, 15251), _mm256_set_epi16(29233, 20476, 26633, 7338, 1495, 14483, 8387, 13130, 32754, 13273, 24223, 24509, 9388, 26364, 21026, 1149), _mm256_set_epi16(19274, 6299, 32683, 3398, 26833, 16830, 2391, 29511, 29233, 20476, 26633, 7338, 1495, 14483, 8387, 13130), _mm256_set_epi16(23902, 29998, 15515, 22406, 11404, 1734, 18742, 24898, 19274, 6299, 32683, 3398, 26833, 16830, 2391, 29511), _mm256_set_epi16(27513, 11637, 8239, 5338, 27752, 8558, 10949, 8271, 23902, 29998, 15515, 22406, 11404, 1734, 18742, 24898), _mm256_set_epi16(9697, 13409, 26552, 22672, 9190, 14438, 12681, 24442, 27513, 11637, 8239, 5338, 27752, 8558, 10949, 8271), _mm256_set_epi16(14626, 22870, 9883, 1352, 27815, 2600, 28964, 18426, 9697, 13409, 26552, 22672, 9190, 14438, 12681, 24442), _mm256_set_epi16(2278, 16787, 31205, 28801, 25983, 30676, 27160, 16897, 14626, 22870, 9883, 1352, 27815, 2600, 28964, 18426)];
    // Layer 5
    for i in 0..8{
        tmp = barrett_mul_26113(a[i+8], twiddle5, twiddle5_r_overq);
        a0[i] = barrett_fake_26113(_mm256_add_epi16(a[i], tmp));
        a0[i+8] = barrett_fake_26113(_mm256_sub_epi16(a[i], tmp));
    }
    // Layer 6
    for i in 0..4{
        tmp = barrett_mul_26113(a[i+4], twiddle6[0], twiddle6_r_overq[0]);
        a1[i] = barrett_fake_26113(_mm256_add_epi16(a0[i], tmp));
        a1[i+4] = barrett_fake_26113(_mm256_sub_epi16(a0[i], tmp));
        tmp = barrett_mul_26113(a[i+12], twiddle6[1], twiddle6_r_overq[1]);
        a1[i+8] = barrett_fake_26113(_mm256_add_epi16(a0[i+8], tmp));
        a1[i+12] = barrett_fake_26113(_mm256_sub_epi16(a0[i+8], tmp));
    }
    // Layer 7
    for i in 0..2{
        for j in 0..4{
            tmp = barrett_mul_26113(a[j*4+i+2], twiddle7[j], twiddle7_r_overq[j]);
            a0[j*4+i] = barrett_fake_26113(_mm256_add_epi16(a1[j*4+i], tmp));
            a0[j*4+i+2] = barrett_fake_26113(_mm256_sub_epi16(a1[j*4+i], tmp));
        }
    }
    // Layer 8
    for i in 0..8{
        tmp = barrett_mul_26113(a[2*i+1], twiddle8[i], twiddle8_r_overq[i]);
        a1[2*i] = barrett_fake_26113(_mm256_add_epi16(a0[2*i], tmp));
        a1[2*i+1] = barrett_fake_26113(_mm256_sub_epi16(a0[2*i], tmp));
    }
    a1
}


pub unsafe fn intt5678_26113(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut tmp: __m256i = _mm256_setzero_si256();
    let twiddle0: [__m256i;8] = [9604, 8171, 8888, 8073, 12760, 2630, 8964, 12314, 8341, 1944, 2138, 7295, 4994, 2380, 10892, 9952];
    let twiddle_r_overq: [[__m256i;8];4];
    // inverse Layer 8
    for i in 0..8{
        a0[2*i] = barrett_fake_26113(_mm256_add_epi16(a[2*i], a[2*i+1]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a[2*i], a[2*i+1]));
        a0[2*i+1] = barrett_mul_26113(tmp, twiddle[3][i], twiddle_r_overq[3][i]);
    }
    // inverse Layer 7
    for i in 0..2{
        for j in 0..4{
            a1[4*j+i] = barrett_fake_26113(_mm256_add_epi16(a0[4*j+i], a0[4*j+i+2]));
            tmp = barrett_fake_26113(_mm256_sub_epi16(a0[4*j+i], a0[4*j+i+2]));
            a1[4*j+i+2] = barrett_mul_26113(tmp, twiddle[2][j*2+i], twiddle_r_overq[2][j*2+i]);
        }
    }
    // inverse Layer 6
    for i in 0..4{
        a0[i] = barrett_fake_26113(_mm256_add_epi16(a1[i], a1[i+4]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a1[i], a1[i+4]));
        a0[i+4] = barrett_mul_26113(tmp, twiddle[1][i], twiddle_r_overq[1][i]);
        a0[i+8] = barrett_fake_26113(_mm256_add_epi16(a1[i+8], a1[i+12]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a1[i+8], a1[i+12]));
        a0[i+12] = barrett_mul_26113(tmp, twiddle[1][i+4], twiddle_r_overq[1][i+4]);
    }
    // inverse Layer 5
    for i in 0..8{
        a1[i] = barrett_fake_26113(_mm256_add_epi16(a1[i], a1[i+8]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a1[i], a1[i+8]));
        a1[i+8] = barrett_mul_26113(tmp, twiddle[0][0], twiddle_r_overq[0][0]);
    }
    a1
}

// inv twiddle factor: [24844, 7249, 21808, 23992, 4719, 20413, 22518, 23068]
//                     [17468, 18908, 20685, 24203] 
//                     [619, 7720]                                     17579]


pub unsafe fn intt1234_26113(a: [__m256i;16]) -> [__m256i;16]{
    let mut a0: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut a1: [__m256i;16] = [_mm256_setzero_si256(); 16];
    let mut tmp: __m256i = _mm256_setzero_si256();
    // inverse Layer 4
    a0[0] = barrett_fake_26113(_mm256_add_epi16(a[0], a[1]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[0], a[1]));
    a0[1] = barrett_mul_24844_26113(tmp);
    a0[2] = barrett_fake_26113(_mm256_add_epi16(a[2], a[3]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[2], a[3]));
    a0[3] = barrett_mul_7249_26113(tmp);
    a0[4] = barrett_fake_26113(_mm256_add_epi16(a[4], a[5]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[4], a[5]));
    a0[5] = barrett_mul_21808_26113(tmp);
    a0[6] = barrett_fake_26113(_mm256_add_epi16(a[6], a[7]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[6], a[7]));
    a0[7] = barrett_mul_23992_26113(tmp);
    a0[8] = barrett_fake_26113(_mm256_add_epi16(a[8], a[9]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[8], a[9]));
    a0[9] = barrett_mul_4719_26113(tmp);
    a0[10] = barrett_fake_26113(_mm256_add_epi16(a[10], a[11]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[10], a[11]));
    a0[11] = barrett_mul_20413_26113(tmp);
    a0[12] = barrett_fake_26113(_mm256_add_epi16(a[12], a[13]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[12], a[13]));
    a0[13] = barrett_mul_22518_26113(tmp);
    a0[14] = barrett_fake_26113(_mm256_add_epi16(a[14], a[15]));
    tmp = barrett_fake_26113(_mm256_sub_epi16(a[14], a[15]));
    a0[15] = barrett_mul_23068_26113(tmp);
    // inverse Layer 3 [17468, 18908, 20685, 24203] 
    for i in 0..2{
        a1[i] = barrett_fake_26113(_mm256_add_epi16(a0[i], a0[i+2]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a0[i], a0[i+2]));
        a1[i+2] = barrett_mul_17468_26113(tmp);
        a1[i+4] = barrett_fake_26113(_mm256_add_epi16(a0[i+4], a0[i+6]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a0[i+4], a0[i+6]));
        a1[i+6] = barrett_mul_18908_26113(tmp);
        a1[i+8] = barrett_fake_26113(_mm256_add_epi16(a0[i+8], a0[i+10]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a0[i+8], a0[i+10]));
        a1[i+10] = barrett_mul_20685_26113(tmp);
        a1[i+12] = barrett_fake_26113(_mm256_add_epi16(a0[i+12], a0[i+14]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a0[i+12], a0[i+14]));
        a1[i+14] = barrett_mul_24203_26113(tmp);
    }
    // inverse Layer 2 [619, 7720]
    for i in 0..4{
        a0[i] = barrett_fake_26113(_mm256_add_epi16(a1[i], a1[i+4]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a1[i], a1[i+4]));
        a0[i+4] = barrett_mul_619_26113(tmp);
        a0[i+8] = barrett_fake_26113(_mm256_add_epi16(a1[i+8], a1[i+12]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a1[i+8], a1[i+12]));
        a0[i+12] = barrett_mul_7720_26113(tmp);
    }
    // inverse Layer 1 [17579]
    for i in 0..8{
        a1[i] = barrett_fake_26113(_mm256_add_epi16(a0[i], a0[i+8]));
        tmp = barrett_fake_26113(_mm256_sub_epi16(a0[i], a0[i+8]));
        a1[i+8] = barrett_mul_17579_26113(tmp);
    }
    a1
}