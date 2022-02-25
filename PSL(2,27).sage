from tqdm import tqdm

_.<x> = QQ[]
Qf.<a> = NumberField(x^12 + 142*x^10 + 7173*x^8 + 157368*x^6 + 1510016*x^4 + 5098688*x^2 + 90352)
Ff.<b> = Qf.subfields()[1][0] 

QQf.<c,d> = Qf.relativize(a^2)
from_QQf, to_QQf = QQf.structure()


l = 3
ll = QQf.base_field().prime_above(l)
lll = QQf.prime_above(ll)

r = lll.absolute_norm()
rprime = ll.absolute_norm()

print("rprime = " + str(rprime))
print("r = " + str(r))

Fl = QQf.base_field().residue_field(ll)
FL = QQf.residue_field(lll)


load("43.5.b.b.qexp.sage")
f = make_data()

print("a2 mod lll = " + str(FL(to_QQf(f[2]))))
print("a7 mod lll = " + str(FL(to_QQf(f[7]))))
print("a2 mod lll in F3^3 = " + str(to_QQf(f[2]) in Fl))

a2 = FL(f[2])
det2 = FL(-2^4)

Ma = 1
Md = a2 - Ma

Mb = 1
Mc = (Ma*Md - det2)/Mb


M = matrix([[Ma, Mb],
            [Mc, Md]])

print("\\bar\\rho_3(Frob_2) = ")
print(M)

i = 1
while M^i != identity_matrix(FL, 2):
    i += 1
print(f"Order of \\bar\\rho_3(Frob_2) = {i}")