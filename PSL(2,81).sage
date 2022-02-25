load("226.2.a.d.qexp.sage")
f = make_data()

Qf = f[5].parent()
FL = Qf.quotient(3)

R.<x> = GF(3)[]

a5 = f[5]
p = a5.minpoly()
print("p_a5 = " + str(p))
print(p.factor())
print("p_a5 % 3 = " + str(R(p)))
print(R(p).factor())
print("a2 mod 3 = " + str(FL(f[2])))