# def barrett_fake(a, f):
#     r_q = round(2^32/f)
#     d = mulh_32(a, r_q)
#     e = mull_32(d, f)
#     return a + f - e

fields = [759207937, 759304193]
j = open("head.txt").read()

for f in fields:
    rq = round(2**32/f)
    print(j)
    print("pub unsafe fn barrett_fake_"+str(f)+"(x: __m512i) -> __m512i{")
    print("    let d = _mm512_mulhi_epu32(x, _mm512_set1_epi32("+str(rq)+"));")
    print("    let e = _mm512_mullo_epi32(d, _mm512_set1_epi32("+str(f)+"));")
    print("    _mm512_sub_epi32(_mm512_add_epi32(x, _mm512_set1_epi32("+str(f)+")), e)\n}\n")