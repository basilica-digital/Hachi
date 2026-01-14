"""
fields = [7681, 10753, 11777, 12289, 13313, 15361]

for f in fields:
    k.<x> = GF(f)
    tmp = (1/f)%(2^16)
    if tmp > 2^15-1:
        tmp -= 2^16
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn montproduct_"+str(f)+"(x: __m256i, y: __m256i) -> __m256i{")
    print("    let lo = _mm256_mullo_epi16(x, y);")
    print("    let hi = _mm256_mulhi_epi16(x, y);")
    print("    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16("+str(tmp)+"));")
    print("    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16("+str(f)+"));")
    print("    let f = _mm256_sub_epi16(hi, e);")
    print("    range_narrow(f, "+str(f)+")")
    print("}\n")
"""