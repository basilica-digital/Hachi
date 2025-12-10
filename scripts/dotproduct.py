field = [7681, 10753, 11777, 12289, 13313, 15361]

for f in field:
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn dot_product_"+str(f)+"(a: [__m256i;16], b: [__m256i;16]) -> [__m256i;16]{")
    print("    let mut c: [__m256i;16] = [_mm256_setzero_si256(); 16]; ")
    print("    for i in 0..16{")
    print("        c[i] = montproduct_"+str(f)+"(a[i], b[i]);")
    print("    }")
    print("    c")
    print("}\n")