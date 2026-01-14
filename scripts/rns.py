"""
f = [7681, 10753, 11777, 12289, 13313, 15361]

def adjust(a, b):
    k.<x> = GF(b)
    tmp = int(1/k(a))
    while tmp > b/2:
        tmp -= b
    while tmp < -b/2:
        tmp += b
    return tmp, round(tmp*2^16/b)

print("for k in 0..16{")
for i in range(1, 6):
    for j in range(i, 6):
        tmp, tmp_r = adjust(f[i-1], f[j])
        print("    v["+str(j)+"][k] = barrett_mul_"+str(f[j])+"(_mm256_sub_epi16(v["+str(j)+"][k], v["+str(i-1)+"][k]), _mm256_set1_epi16("+str(tmp)+"), _mm256_set1_epi16("+str(tmp_r)+"));")
print("}\n")
"""