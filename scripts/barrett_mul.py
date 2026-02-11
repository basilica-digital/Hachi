# def barrett_mul(a, b, f):
#     neg_f = -f & (2^32-1)
#     br = floor(b*2^32/f)
#     t = mulh_32(a, br)
#     d = mull_32(a, b)
#     g = mull_32(t, neg_f)
#     return (d + g)&(2^32-1)

fields = [759207937, 759304193]
j = open("head.txt").read()

for f in fields:
    rq = round(2**32/f)
    print(j)
    print("pub unsafe fn barrett_mul"+str(f)+"(x: __m512u) -> __m512u{")
    print("    let d = _mm512_mulhi_epu32(x, _mm512_set1_epu32("+str(rq)+"));")
    print("    let e = _mm512_mullo_epi32(d, _mm512_set1_epu32("+str(f)+"));")
    print("    _mm256_sub_epi16(x, e);\n}\n")