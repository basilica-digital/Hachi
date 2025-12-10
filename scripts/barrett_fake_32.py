f = [7681, 10753, 11777, 12289, 13313, 15361]
b = [5585133, 1290167, 559168, 541132, 453390, 399420]

for i in range(len(f)):
    print("#[target_feature(enable = \"avx2\")]")
    print("pub unsafe fn barrett_fake_"+str(f[i])+"_32(x: __m256i, y: __m256i) -> __m256i{")
    print("    let m1 = _mm256_set1_epi32("+str(b[i])+");")
    print("    let m2 = _mm256_set1_epi32("+str(f[i])+");")
    print("    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32("+str(b[i])+"));")
    print("    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32("+str(f[i])+"));")
    print("    let f = _mm256_sub_epi16(x, e);")
    print("    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32("+str(b[i])+"));")
    print("    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32("+str(f[i])+"));")
    print("    let f1 = _mm256_sub_epi16(y, e1);")
    print("    let ans = pack_32_16(f, f1);\n}\n")
    