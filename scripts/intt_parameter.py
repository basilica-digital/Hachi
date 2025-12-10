"""
fields = [7681, 10753, 11777, 12289, 13313, 15361]

def root_of_unity_to_mm256(z, f):
    k.<x> = GF(f)
    a = []
    for i in z:
        a.append(int(1/(k(i).sqrt())))
    b = [[] for i in range(len(a)//16)]
    for i in range(len(b)):
        for j in range(16):
            b[i].append(a[i*16+j])
            
    s = ""
    for i in b:
        s += "_mm256_set_epi16("
        for j in range(len(i)-1, 0, -1):
            s += str(i[j]) + ", "
        s = s[:-2] + "), _mm256_set_epi16("
        for j in range(len(i)-1, 0, -1):
            s += str(round(i[j]*2^16/f)) + ", "
        s = s[:-2] + "), "
    return s[:-2]
    

for f in fields:
    k.<x> = GF(f)
    tmp = [[]for i in range(8)]
    tmp[0].append(int(k(-1).sqrt()))
    tmp[0].append(-tmp[0][0])
    print("    let l1:[__m256i;2] = [_mm256_set1_epi16("+str(int(k(-1).sqrt()))+"), _mm256_set1_epi16("+str(round(int(k(-1).sqrt())*2^16/f))+")];")
    for i in range(1, 4):
        s = "    let l" + str(i+1) + ":[__m256i;"+str(2^(i+1))+"] = ["
        for j in range(len(tmp[i-1])):
            t = int(1/(k(tmp[i-1][j]).sqrt()))
            s += "_mm256_set1_epi16("+str(t)+"), "
            s += "_mm256_set1_epi16("+str(round(t*2^16/f))+"), "
            tmp[i].append(t)
            tmp[i].append(-t)
        print(s[:-2] + "];")
    print("\n\n")
    
    for i in range(4, 8):
        print("    let l" + str(i+1) + ":[__m256i;"+str(2^(i-3))+"] = [" + root_of_unity_to_mm256(tmp[i-1], f) + "];")
        for j in range(len(tmp[i-1])):
            t = int(k(tmp[i-1][j]).sqrt())
            tmp[i].append(t)
            tmp[i].append(-t)
    print("\n\n")
"""