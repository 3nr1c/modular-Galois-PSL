_.<x> = QQ[]
Qf.<a> = NumberField(x^12 + 7208 * x^10 + 19859688 * x^8 + 26566749360 * x^6 + 17884354852944 * x^4 + 5570285336959680 * x^2 + 590986232936064000)
Ff.<b> = Qf.subfields()[1][0] 

print(Ff.degree())

l = 7

QQf.<c,d> = Qf.relativize(a^2)
ll = QQf.base_field().primes_above(l)[2]
lll = QQf.prime_above(ll)

from_QQf, to_QQf = QQf.structure()

print(ll.absolute_norm().factor())
print(lll.absolute_norm().factor())

Fl = QQf.base_field().residue_field(ll)
FL = QQf.residue_field(lll)


load("31.7.b.c.qexp.sage")
f = make_data()

print("a3 mod lll in F7^3 = " + str(to_QQf(f[3]) in Fl))
print("a3 mod lll = " + str(FL(to_QQf(f[3]))))
print("a5 mod lll = " + str(FL(to_QQf(f[5]))))
print("a7 mod lll = " + str(FL(to_QQf(f[7]))))


a3 = FL(to_QQf(f[3]))
det3 = FL(- 3^6)

a = 1
d = a3 - 1

b = 1
c = (a*d - det3)/b

M = matrix([[a, b],
            [c, d]])

print("\\bar\\rho_7(Frob_3) = ")
print(M)

for i in range(1, 1000):
    if M^i == identity_matrix(FL, 2):
        print(f"Order of \\bar\\rho_7(Frob_3) = {i}")
        break