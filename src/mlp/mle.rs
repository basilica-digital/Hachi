use crate::math::fields::*;
use crate::utils::ds::*;
//move
use crate::Q;

pub type MBundle = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);
pub type WBundle = (AlignedU32Vec, AlignedU32Vec, AlignedU32Vec, AlignedU32Vec, AlignedU32Vec, AlignedI8Vec);


fn mle_M(i: u32, u: u32, m: &MBundle) -> [u32; 4]{
    let zero = [0u32, 0u32, 0u32, 0u32];
    let u = u as usize;
    let i = i as usize;
    // D
    if i == 0{
        if u < (1<<13){
            return m.0[4*u..4*(u+1)].try_into().unwrap();
        }else if u < (2<<13) && u >= (1<<13){
            let idx = u - (1<<13);
            return m.13[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (3<<13) && u >= (2<<13){
            let idx = u - (2<<13);
            return m.14[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (4<<13) && u >= (3<<13){
            let idx = u - (3<<13);
            return m.15[4*idx..4*(idx+1)].try_into().unwrap();
        }else {
            return zero;
        }
    }
    // B
    if i == 1{
        if u < (5<<13) && u >= (4<<13){
            let idx = u - (4<<13);
            return m.1[4*idx..4*(idx+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    // bG
    else if i == 2{
        if u < 1<<13{
            return m.7[4*u..4*(u+1)].try_into().unwrap();
        }else if u < (2<<13) && u >= (1<<13){
            let idx = u - (1<<13);
            return extdbl(m.10[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (3<<13) && u >= (2<<13){
            let idx = u - (2<<13);
            return extdbl(m.9[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (4<<13) && u >= (3<<13){
            let idx = u - (3<<13);
            return extdbl(m.8[4*idx..4*(idx+1)].try_into().unwrap());
        }else{
            return zero;
        }
    }
    else if i == 3{
        if u < 1<<13{
            return m.8[4*u..4*(u+1)].try_into().unwrap();
        }else if u < (2<<13) && u >= (1<<13){
            let idx = u - (1<<13);
            return m.7[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (3<<13) && u >= (2<<13){
            let idx = u - (2<<13);
            return extdbl(m.10[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (4<<13) && u >= (3<<13){
            let idx = u - (3<<13);
            return extdbl(m.9[4*idx..4*(idx+1)].try_into().unwrap());
        }else{
            return zero;
        }
    }
    else if i == 4{
        if u < 1<<13{
            return m.9[4*u..4*(u+1)].try_into().unwrap();
        }else if u < (2<<13) && u >= (1<<13){
            let idx = u - (1<<13);
            return m.8[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (3<<13) && u >= (2<<13){
            let idx = u - (2<<13);
            return m.7[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (4<<13) && u >= (3<<13){
            let idx = u - (3<<13);
            return extdbl(m.10[4*idx..4*(idx+1)].try_into().unwrap());
        }else{
            return zero;
        }
    }
    else if i == 5{
        if u < 1<<13{
            return m.10[4*u..4*(u+1)].try_into().unwrap();
        }else if u < (2<<13) && u >= (1<<13){
            let idx = u - (1<<13);
            return m.9[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (3<<13) && u >= (2<<13){
            let idx = u - (2<<13);
            return m.8[4*idx..4*(idx+1)].try_into().unwrap();
        }else if u < (4<<13) && u >= (3<<13){
            let idx = u - (3<<13);
            return m.7[4*idx..4*(idx+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    // cG, aG
    else if i == 6{
        if u >= 5<<13 && u < 9<<13 {
            let idx = u - (5<<13);
            return extneg(m.3[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (1<<13){
            return m.11[4*u..4*(u+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    else if i == 7{
        if u >= 5<<13 && u < 9<<13 {
            let idx = u - (5<<13);
            return extneg(m.4[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (2<<13) && u >= (1<<13){
            let idx = u - (1<<13);
            return m.11[4*idx..4*(idx+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    else if i == 8{
        if u >= 5<<13 && u < 9<<13 {
            let idx = u - (5<<13);
            return extneg(m.5[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (3<<13) && u >= (2<<13){
            let idx = u - (2<<13);
            return m.11[4*idx..4*(idx+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    else if i == 9{
        if u >= 5<<13 && u < 9<<13 {
            let idx = u - (5<<13);
            return extneg(m.6[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (4<<13) && u >= (3<<13){
            let idx = u - (3<<13);
            return m.11[4*idx..4*(idx+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    // cG, A
    else if i == 10{
        if u >= 5<<13 && u < 9<<13 {
            let idx = u - (5<<13);
            return extneg(m.2[4*idx..4*(idx+1)].try_into().unwrap());
        }else if u < (5<<13) && u >= (4<<13){
            let idx = u - (4<<13);
            return m.11[4*idx..4*(idx+1)].try_into().unwrap();
        }else{
            return zero;
        }
    }
    return zero;
}

pub fn mle_m(i: u32, u: u32, m: &MBundle) -> [u32; 4]{
    // inside M
    if (i < 11) && (u < (9<<13)){
        mle_M(i, u, m)
    }
    // -(alpha^1024 + 1)
    else if (u >= 9<<13) && (u < (9<<13)+11) && (i as usize == (u - (9<<13)) as usize) { 
        m.12.clone().try_into().unwrap()
    }
    // 0
    else{ 
        [0u32, 0u32, 0u32, 0u32]
    }
}

fn mle_W(u: u32, l: u32, w: &WBundle) -> u32{
    if u < (1<<13) {
        w.0[u as usize *1024 + l as usize]
    }else if (u < (2<<13) && u >= (1<<13)){
        w.1[(u as usize -(1<<13))*1024+ l as usize]
    }else if (u < (3<<13) && u >= (2<<13)){
        w.2[(u as usize -(2<<13))*1024+ l as usize]
    }else if (u < (4<<13) && u >= (3<<13)){
        w.3[(u as usize -(3<<13))*1024+ l as usize]
    }else if (u < (5<<13) && u >= (4<<13)){
        w.4[(u as usize -(4<<13))*1024+ l as usize]
    }else if (u < (9<<13) && u >= (5<<13)){
        let val = w.5[(u as usize - (5 << 13)) * 1024 + l as usize]; 
        if val < 0 {
            return 4294967197u32 - ((-val) as u32);
        } else {
            return val as u32;
        }
    } else {
        0
    }
}

pub fn mle_w(u: u32, l: u32, w: &WBundle, r: &[u32]) -> u32{
    // inside W
    if u < (9<<13) {
        mle_W(u, l, w)
    }
    // r
    else if u < (9<<13) + 11 { 
        r[(u as usize -(9<<13))*1024 + l as usize]
    } else {
        0
    }
}

pub fn mle_alpha(y: u32, alpha_vec: &[u32]) -> u32{
    alpha_vec[y as usize]
}