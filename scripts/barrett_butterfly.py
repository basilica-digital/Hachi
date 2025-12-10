f = [7681, 10753, 11777, 12289, 13313, 15361]

for i in f:
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn barrett_butterfly_"+str(i)+"(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{")
    print("    let t = _mm256_mulhi_epi16(b, cr);")
    print("    let d = _mm256_mullo_epi16(b, c);")
    print("    let dp = _mm256_add_epi16(a, d);")
    print("    let ds = _mm256_sub_epi16(a, d);")
    print("    let e = _mm256_sub_epi16(t, _mm256_set1_epi16("+str(i)+"));")
    print("    let ansp = _mm256_sub_epi16(dp, e);")
    print("    let anss = _mm256_add_epi16(dp, e);")
    print("    [ansp, anss]")
    print("}\n")

