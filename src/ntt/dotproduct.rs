#[target_feature(enable = "avx2")]
pub unsafe fn dot_product_7681(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{
    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16]; 
    for i in 0..16{
        c[i] = montproduct_7681(a[i], b[i]);
    }
    c
}

#[target_feature(enable = "avx2")]
pub unsafe fn dot_product_10753(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{
    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16]; 
    for i in 0..16{
        c[i] = montproduct_10753(a[i], b[i]);
    }
    c
}

#[target_feature(enable = "avx2")]
pub unsafe fn dot_product_11777(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{
    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16];
    for i in 0..16{
        c[i] = montproduct_11777(a[i], b[i]);
    }
    c
}

#[target_feature(enable = "avx2")]
pub unsafe fn dot_product_12289(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{
    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16];
    for i in 0..16{
        c[i] = montproduct_12289(a[i], b[i]);
    }
    c
}

#[target_feature(enable = "avx2")]
pub unsafe fn dot_product_13313(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{
    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16];
    for i in 0..16{
        c[i] = montproduct_13313(a[i], b[i]);
    }
    c
}

#[target_feature(enable = "avx2")]
pub unsafe fn dot_product_15361(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{
    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16];
    for i in 0..16{
        c[i] = montproduct_15361(a[i], b[i]);
    }
    c
}