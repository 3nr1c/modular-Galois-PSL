from tqdm import tqdm

_.<x> = QQ[]
Qf.<a> = NumberField(x^6 + 398*x^4 + 49236*x^2 + 1934136)
Ff.<b> = Qf.subfields()[1][0] 

QQf.<c,d> = Qf.relativize(a^2)
from_QQf, to_QQf = QQf.structure()

l = 5
ll = QQf.base_field().prime_above(5)
lll = QQf.prime_above(ll)

r = lll.absolute_norm()
rprime = ll.absolute_norm()

print("rprime = " + str(rprime))
print("r = " + str(r))

Fl = QQf.base_field().residue_field(ll)
FL = QQf.residue_field(lll)


load("31.5.b.b.qexp.sage")
f = make_data()

print("a2 mod lll = " + str(FL(to_QQf(f[2]))))
print("a3 mod lll = " + str(FL(to_QQf(f[3]))))
print("a5 mod lll = " + str(FL(to_QQf(f[5]))))
print()

a3 = FL(to_QQf(f[3]))
det3 = FL(-3^4)

a = 1
d = a3 - a

b = 1
c = (a*d - det3)/b

a = 1
d = a3 - 1
b = 1
c = a3

M = matrix([[a, b],
            [c, d]])

print("\\bar\\rho_5(Frob_3) = ")
print(M)
print()

for i in range(1, 976500):
    if M^i == identity_matrix(FL, 2):
        print(f"Order of \\bar\\rho_5(Frob_3) = {i}")
        break

