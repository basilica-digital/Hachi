"""f = [7681, 10753, 11777, 12289, 13313, 15361]

print("for k in 0..16{")
for i in range(1, 6):
    for j in range(i, 6):
        k.<x> = GF(f[j])
        tmp = int(1/k(f[i-1]))
        tmp_r = round(tmp*2^16/f[j])
        print("    v["+str(j)+"][k] = barrett_mul_"+str(f[j])+"(_mm256_sub_epi16(v["+str(j)+"][k], v["+str(i-1)+"]), _mm256_set1_epi16("+str(tmp)+"));")
print("}\n")"""