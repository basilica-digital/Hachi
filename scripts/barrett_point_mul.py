f = [7681, 10753, 11777, 12289, 13313, 15361]
b = [57857, 54785, 53761, 53249, 52225, 50177]

for i in range(len(f)):
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn montproduct_"+str(f[i])+"(x: __m256i, y: __m256i) -> __m256i{")
    print("    let lo = _mm256_mullo_epi16(x, y);")
    print("    let hi = _mm256_mulhi_epi16(x, y);")
    print("    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16("+str(b[i])+"));")
    print("    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16("+str(f[i])+"));")
    print("    let f = _mm256_sub_epi16(hi, e);\n}\n")