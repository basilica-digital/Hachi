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
primes = [7681, 10753, 11777, 12289, 13313, 15361]
N = 2^32 - 99

prodx = [[0, 7681, 13057, 8679, 4808, 0], [0, 7681, 13057, 8679, 4808, 0], [0, 7681, 13057, 8679, 4808, 0], [0, 7681, 11713, 8679, 4808, 0], [0, 7681, 13254, 10746, 4808, 0], [0, 7681, 13057, 8679, 4808, 0]]
constants = [5380, 11074, 512, 9776, 14238]
# v[0] = residues[0]

# 現在是做一個 vector
for i in range(1, 6):
    # subtract_val = v[0]
    print("\tlet mut subtract_val = a;")
    for j in range(1, i):
        # subtract_val = (subtract_val + v[j] * prodx[i][j]) % primes[i]
        print("subtract_val = barrett_mla_"+ str(primes[i]) +"(subtract_val, a[j], _mm256_set1_epi16("+ str(prodx[i][j]) +"), _mm256_set1_epi16("+ str(round(prodx[i][j]*2^16/primes[i])) +"));")
    # v[i] = (((residues[i] - subtract_val)) * constants[i-1]) % primes[i]
    print("v[i] = barrett_mul_"+ str(primes[i]) +"(_mm256_sub_epi16(a[i], subtract_val), "+ str(constants[i-1]) +");")