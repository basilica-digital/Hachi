"""
def to_int16(x):
    u = Integer(x) & 0xFFFF 
    if u >= 0x8000:
        return u - 0x10000
    return u

f = [769, 3329, 7681, 7937, 9473, 10753]
for i in f:
    for j in range(10000):
        r = ZZ.random_element(-(2^31), 2^32-1)
        m1 = round((m*(2^16))/(i))
        t = (r*m1) >> 16
        d = to_int16(r*m)
        dp = to_int16(a+d)
        ds = to_int16(a-d)
        ansp = to_int16(dp - to_int16(t*i))
        anss = to_int16(ds + to_int16(t*i))
        if ansp%i != (a + r*m)%i or ansp - (a + (r*m)%i) < -3*i or ansp - (a + (r*m)%i) > 3*i:
            print(ansp, a + (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
        if anss%i != (a - r*m)%i or anss - (a - (r*m)%i) < -3*i or anss - (a - (r*m)%i) > 3*i:
            print(anss, a - (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
"""