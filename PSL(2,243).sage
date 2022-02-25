_.<x> = QQ[]
Qf.<a> = NumberField(x^10 + 32 * x^8 + 357 * x^6 + 1725 * x^4 + 3366 * x^2 + 1519)
Ff.<b> = Qf.subfields()[1][0] 

print(Ff.degree())

l = 3

QQf.<c,d> = Qf.relativize(a^2)
ll = QQf.base_field().prime_above(l)
lll = QQf.prime_above(ll)


from_QQf, to_QQf = QQf.structure()

print(ll.absolute_norm().factor())
print(lll.absolute_norm().factor())

Fl = QQf.base_field().residue_field(ll)
FL = QQf.residue_field(lll)

load("67.3.b.b.qexp.sage")
f = make_data()

print("a2 mod lll = " + str(FL(to_QQf(f[2]))))
print("a3 mod lll = " + str(FL(to_QQf(f[3]))))
print("a7 mod lll = " + str(FL(to_QQf(f[7]))))


a2 = FL(to_QQf(f[2]))
det2 = FL(-2^2)

a = 1
d = a2 - 1

b = 1
c = (a*d - det2)/b

M = matrix([[a, b],
            [c, d]])

print("\\bar\\rho_3(Frob_2) = ")
print(M)

for i in range(1, 1000):
    if M^i == identity_matrix(FL, 2):
        print(f"Order of \\bar\\rho_3(Frob_2) = {i}")
        break