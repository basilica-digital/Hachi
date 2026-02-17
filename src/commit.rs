#[inline(never)]
pub unsafe fn commit(u: &mut [u32], s: &mut [u32], acap1: &mut [u32], acap2: &mut [u32]){
    let n = 1 << 10;
    let moduli = [759207937, 759304193];
    let plans: Vec<Plan> = moduli.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();
    let iterations = 1 << 10;
    let height = 1 << 13;

    let mut t = vec![0u32; iterations*n];
    let mut gt = vec![0u32; height*n];
    // let mut capa1 = vec![0u32; height*n];
    // let mut capa2 = vec![0u32; height*n];

    let gt_ptr = gt.as_mut_ptr() as *mut __m512i;
    let t_ptr = t.as_mut_ptr() as *mut __m512i;
    let u_ptr = u.as_mut_ptr() as *mut __m512i;
    let capa1_ptr = acap1.as_mut_ptr() as *mut __m512i;
    let capa2_ptr = acap2.as_mut_ptr() as *mut __m512i;

    // // RNS and ntt on A
    // for i in 0..height{
    //     for j in (0..1024).step_by(16) {
    //         rns_decompose_16_elements(&acap[i*1024+j..i*1024+(j+16)], &mut capa1[i*1024+j..i*1024+(j+16)], &mut capa2[i*1024+j..i*1024+(j+16)], j);
    //     }
    //     let mut a1 = &mut capa1[i*1024..(i+1)*1024];
    //     let mut a2 = &mut capa2[i*1024..(i+1)*1024];
    // }

    // As
    for i in 0..iterations{
        let mut val1 = [_mm512_set1_epi32(0); 64];
        let mut val2 = [_mm512_set1_epi32(0); 64];
        for j in 0..height{
            let mut s1:[u32;1024] = (&s[j*1024+i..(j+1)*1024+i]).try_into().expect("Slice length mismatch");
            let mut s2:[u32;1024] = (&s[j*1024+i..(j+1)*1024+i]).try_into().expect("Slice length mismatch");
            plans[0].fwd(&mut s1);
            plans[1].fwd(&mut s2);
            let s1_ptr = s1.as_mut_ptr() as *mut __m512i;
            let s2_ptr = s2.as_mut_ptr() as *mut __m512i;
            for k in 0..64{
                let s1_vec = _mm512_loadu_si512(s1_ptr.add(k));
                let s2_vec = _mm512_loadu_si512(s2_ptr.add(k));
                let a1_vec = _mm512_loadu_si512(capa1_ptr.add(k));
                let a2_vec = _mm512_loadu_si512(capa2_ptr.add(k));
                val1[k] = barrett_fake_759207937(_mm512_add_epi32(val1[k], montproduct_759207937(s1_vec, a1_vec)));
                val2[k] = barrett_fake_759304193(_mm512_add_epi32(val2[k], montproduct_759304193(s2_vec, a2_vec)));  
            }
        }
        plans[0].inv(slice::from_raw_parts_mut(val1.as_mut_ptr() as *mut u32, 1024));
        plans[1].inv(slice::from_raw_parts_mut(val2.as_mut_ptr() as *mut u32, 1024));
        irns(&val1, &val2, t_ptr, i);
    }

    // decompose t
    let mask15 = _mm512_set1_epi32(15);
    for i in 0..iterations{
        for j in 0..64{
            let mut val = _mm512_loadu_si512(t_ptr.add(i*64+j));
            for k in 0..8{
                _mm512_storeu_si512((gt_ptr.add(i*512+k*64+j)) as *mut _ , _mm512_and_epi32(mask15, val));
                val = _mm512_srli_epi32(val, 4);
            }
        }
    }
    let mut v1 = [_mm512_set1_epi32(0); 64];
    let mut v2 = [_mm512_set1_epi32(0); 64];
    for i in 0..height{
        let mut gt1:[u32;1024] = (&gt[i*1024..(i+1)*1024]).try_into().expect("Slice length mismatch");
        let mut gt2:[u32;1024] = (&gt[i*1024..(i+1)*1024]).try_into().expect("Slice length mismatch");
        plans[0].fwd(&mut gt1);
        plans[1].fwd(&mut gt2);
        let gt1_ptr = gt1.as_mut_ptr() as *mut __m512i;
        let gt2_ptr = gt2.as_mut_ptr() as *mut __m512i;
        for j in 0..64{
            let gt1_vec = _mm512_loadu_si512(gt1_ptr.add(j));
            let gt2_vec = _mm512_loadu_si512(gt2_ptr.add(j));
            let b1_vec = _mm512_loadu_si512(capa1_ptr.add(j));
            let b2_vec = _mm512_loadu_si512(capa2_ptr.add(j));
            v1[j] = barrett_fake_759207937(_mm512_add_epi32(v1[j], montproduct_759207937(gt1_vec, b1_vec)));
            v2[j] = barrett_fake_759304193(_mm512_add_epi32(v2[j], montproduct_759304193(gt2_vec, b2_vec)));  
        }
    }
    plans[0].inv(slice::from_raw_parts_mut(v1.as_mut_ptr() as *mut u32, 1024));
    plans[1].inv(slice::from_raw_parts_mut(v2.as_mut_ptr() as *mut u32, 1024));
    irns(&v1, &v2, u_ptr, 0);
}