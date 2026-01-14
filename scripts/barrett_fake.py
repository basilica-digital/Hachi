f = [7681, 10753, 11777, 12289, 13313, 15361]
b = [round((2**16)/i) for i in f]

for i in range(len(f)):
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn barrett_fake_"+str(f[i])+"(x: __m256i) -> __m256i{")
    print("    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16("+str(b[i])+"));")
    print("    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16("+str(f[i])+"));")
    print("    let f = _mm256_sub_epi16(x, e);\n    f\n}\n")