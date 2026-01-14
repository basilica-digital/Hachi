f = [7681, 10753, 11777, 12289, 13313, 15361]

for i in f:
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn barrett_mul_"+str(i)+"(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{")
    print("    let b = range_narrow(bb, "+str(i)+");")
    print("    let t = _mm256_mulhi_epi16(b, ar_overq);")
    print("    let d = _mm256_mullo_epi16(b, a);")
    print("    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16("+str(i)+"));")
    print("    let g = _mm256_sub_epi16(d, f);")
    print("    range_narrow(g, "+str(i)+")\n}\n")