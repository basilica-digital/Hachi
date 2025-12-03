"""
fields = [257, 3329, 7681, 7937, 9473, 10753, 26113]

for f in fields:
    k.<x> = GF(f)
    tbl = [[]for i in range(4)]
    tmp = [[]for i in range(4)]
    tmp[0].append(int(k(-1).sqrt()))
    tmp[0].append(-tmp[0][0])
    tbl[0].append(int(k(-1).sqrt()))
    tbl[0].append(round(tbl[0][0]*2^16/f))
    print("    let l0:[__m256i;2] = [_mm256_set1_epi16("+str(tbl[0][0])+"), _mm256_set1_epi16("+str(tbl[0][1])+")];")
    for i in range(1, 4):
        s = "    let l" + str(i) + ":[__m256i;"+str(2^(i+1))+"] = ["
        for j in range(len(tmp[i-1])):
            t = int(k(tmp[i-1][j]).sqrt())
            tbl[i].append(t)
            s += "_mm256_set1_epi16("+str(t)+"), "
            s += "_mm256_set1_epi16("+str(round(t*2^16/f))+"), "
            tbl[i].append(round(t*2^16/f))
            tmp[i].append(t)
            tmp[i].append(-t)
        print(s[:-2] + "];")
    print("\n\n")
"""