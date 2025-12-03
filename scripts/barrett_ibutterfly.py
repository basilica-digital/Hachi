f = [257, 3329, 7681, 7937, 9473, 10753]

for i in f:
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn barrett_ibutterfly_"+str(i)+"(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{")
    print("    let ans0 = barrett_fake_"+str(i)+"(_mm256_add_epi16(a, b));")
    print("    let tmp = _mm256_sub_epi16(a, b);")
    print("    let ans1 = barrett_mul(tmp, c, cr);")
    print("    [ans0, ans1]")
    print("}\n")

