"""
l = []
for i in range(1, 256):
    a = 2^8*(i)+1
    if a > 2^8 and is_prime(a):
        l.append(a)
print(l)
t = [769, 3329, 7681, 7937, 9473, 10753]
z = 1
for k in t:
    z = z*k
print(z - 256*((2^32-99)^2))
"""