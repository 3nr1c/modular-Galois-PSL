
# q-expansion of newform 226.2.a.d, downloaded from the LMFDB on 07 February 2022.

# We generate the q-expansion using the Hecke eigenvalues a_p at the primes.
# Each a_p is given as a linear combination
# of the following basis for the coefficient ring.
# To create the q-expansion as a power series, type "qexp = make_data()"

def make_data():
    from sage.all import prod, floor, prime_powers, gcd, QQ, primes_first_n, next_prime, RR

    def discrete_log(elts, gens, mod):
        # algorithm 2.2, page 16 of https://arxiv.org/abs/0903.2785
        def table_gens(gens, mod):
            T = [1]
            n = len(gens)
            r = [None]*n
            s = [None]*n
            for i in range(n):
                beta = gens[i]
                r[i] = 1
                N = len(T)
                while beta not in T:
                    for Tj in T[:N]:
                        T.append((beta*Tj) % mod)
                    beta = (beta*gens[i]) % mod
                    r[i] += 1
                s[i] = T.index(beta)
            return T, r, s
        T, r, s = table_gens(gens, mod)
        n = len(gens)
        N = [ prod(r[:j]) for j in range(n) ]
        Z = lambda s: [ (floor(s/N[j]) % r[j]) for j in range(n)]
        return [Z(T.index(elt % mod)) for elt in elts]
    def extend_multiplicatively(an):
        for pp in prime_powers(len(an)-1):
            for k in range(1, (len(an) - 1)//pp + 1):
                if gcd(k, pp) == 1:
                    an[pp*k] = an[pp]*an[k]
    from sage.all import PolynomialRing, NumberField, ZZ
    R = PolynomialRing(QQ, "x")
    f = R(poly_data)
    K = NumberField(f, "a")
    betas = [K([c/ZZ(den) for c in num]) for num, den in basis_data]
    convert_elt_to_field = lambda elt: sum(c*beta for c, beta in zip(elt, betas))
    # convert aps to K elements
    primes = primes_first_n(len(aps_data))
    good_primes = [p for p in primes if not p.divides(level)]
    aps = map(convert_elt_to_field, aps_data)
    if not hecke_ring_character_values:
        # trivial character
        char_values = dict(zip(good_primes, [1]*len(good_primes)))
    else:
        gens = [elt[0] for elt in hecke_ring_character_values]
        gens_values = [convert_elt_to_field(elt[1]) for elt in hecke_ring_character_values]
        char_values = dict([(
            p,prod(g**k for g, k in zip(gens_values, elt)))
            for p, elt in zip(good_primes, discrete_log(good_primes, gens, level))
            ])
    an_list_bound = next_prime(primes[-1])
    an = [0]*an_list_bound
    an[1] = 1
    
    from sage.all import PowerSeriesRing
    PS = PowerSeriesRing(K, "q")
    for p, ap in zip(primes, aps):
        if p.divides(level):
            euler_factor = [1, -ap]
        else:
            euler_factor = [1, -ap, p**(weight - 1) * char_values[p]]
        k = RR(an_list_bound).log(p).floor() + 1
        foo = (1/PS(euler_factor)).padded_list(k)
        for i in range(1, k):
            an[p**i] = foo[i]
    extend_multiplicatively(an)
    return PS(an)
level  =  226
weight  =  2
poly_data  =  [5, 0, -5, 0, 1]

# The entries in the following list give a basis for the
# coefficient ring in terms of a root of the defining polynomial above.
# Each line consists of the coefficients of the numerator, and a denominator.
basis_data   =  [\
[[1, 0, 0, 0], 1],
[[-3, 1, 1, 0], 1],
[[3, 1, -1, 0], 1],
[[3, -3, -1, 1], 1]]

hecke_ring_character_values  =  None
aps_data  =  [\
[1, 0, 0, 0],
[1, 1, -1, 1],
[1, -1, 0, -1],
[-2, -1, 1, 0],
[0, 1, 1, 0],
[2, 0, 0, -2],
[-2, -2, 2, 0],
[-3, -1, 1, 1],
[0, 0, -3, 0],
[-1, -1, 0, 1],
[2, 3, -1, 0],
[-3, -3, 0, -1],
[2, 3, 1, 2],
[-1, 1, -3, 3],
[0, 2, -1, 0],
[0, 2, 0, 4],
[7, 1, -1, -1],
[4, 2, 4, -2],
[-5, -1, 1, -1],
[8, 2, -3, 0],
[0, 2, -4, 6],
[-2, -6, 3, -6],
[2, -1, 3, -2],
[14, 0, 0, 0],
[6, 4, -6, 6],
[3, -3, -4, 1],
[-8, -2, -3, 4],
[1, -1, 3, -3],
[0, 4, -2, 4],
[-1, 0, 0, 0],
[-6, 1, 3, -4],
[2, -4, -4, 2],
[-2, 4, -4, 8],
[-8, -7, 5, -4],
[4, -10, 4, -6],
[-8, 0, -1, 0],
[2, 8, -4, 6],
[-2, -4, 4, -2],
[2, -6, 5, -10],
[0, 0, 2, 0],
[-1, 3, 1, 7],
[-3, -3, 4, -5],
[-4, 0, 7, -4],
[-8, -6, 4, -2],
[9, 5, -8, -1],
[-2, 8, -5, 6],
[-4, 2, -6, 4],
[0, 2, -1, -4],
[-4, 4, 0, 0],
[1, 1, -4, -1],
[-2, -1, 1, -6],
[6, -5, 7, -8],
[-6, -5, 7, -4],
[-6, -8, 8, 2],
[-10, -8, 2, 2],
[2, 6, -5, 2],
[17, 5, 0, 3],
[-2, 2, -1, -6],
[8, 6, 0, -6],
[-6, 2, 10, -4],
[6, -5, 3, -10],
[-5, 1, 0, 1],
[4, 0, 0, 0],
[10, -1, -7, 0],
[18, 5, 1, 4],
[-2, 6, 2, 2],
[14, -1, -1, 2],
[6, 7, 9, -6],
[-16, -6, 6, 0],
[-7, -1, 4, 11],
[-18, -7, 7, -6],
[12, 4, -11, 8],
[-12, -10, 10, 4],
[13, -1, 0, 3],
[-7, -5, 9, 5],
[-6, -1, -5, 8],
[-2, 14, -6, 2],
[1, 5, 4, -1],
[2, -3, -7, -4],
[-10, 2, -6, 8],
[23, 3, 5, -9],
[-8, -6, 4, -2],
[-18, -8, 15, -14],
[-8, 4, 2, -6],
[-4, -6, 12, 4],
[-16, -1, 3, -8],
[6, -4, 4, -4],
[-4, -10, 12, -6],
[-20, -6, 4, 2],
[-20, -2, 4, -4],
[2, -5, 11, -14],
[2, 0, -9, 2],
[-4, -2, 3, -4],
[7, 11, -11, -1],
[3, 5, -11, 15],
[-20, -6, 8, -4],
[2, 4, -12, 8],
[10, 0, -2, 6],
[31, 1, 1, -1],
[-1, 5, -4, 17],
[-4, 11, 3, 4],
[4, 2, -12, 0],
[-4, -12, 0, -8],
[6, -8, 6, -2],
[1, -7, 3, 1],
[-2, -2, 2, 8],
[-18, -4, 0, 6],
[10, 8, -10, 6],
[22, 16, -19, 14],
[30, 0, 6, -2],
[-16, 8, -7, 8],
[-3, 1, 8, 7],
[6, 15, -15, 6],
[3, -11, 9, -9],
[16, -2, -3, -4],
[-4, 10, -8, 6],
[5, 11, -5, -3],
[4, -10, 8, -12],
[-4, 6, -4, -2],
[-9, 5, 5, 3],
[1, -7, 4, -9],
[-20, -8, 10, -10],
[22, 16, -16, 16],
[9, 11, 3, -7],
[2, 10, -10, 14],
[-5, 5, 0, -3],
[22, -8, 0, -2],
[26, -1, 7, -8],
[-14, 1, -13, -8],
[13, -1, -4, -1],
[-20, -1, -5, 16],
[-18, -12, 11, 2],
[-2, 8, 5, 2],
[35, 3, -4, -7],
[22, 3, 3, -16],
[-10, -1, -7, -10],
[12, 2, -20, 14],
[14, 19, -5, -2],
[3, -19, 8, -7],
[-6, -15, 11, 2],
[19, 9, -7, -1],
[6, -4, -4, 8],
[-16, 14, -6, 16],
[0, -9, 7, 0],
[-13, -3, 20, -3],
[2, -8, 1, -22],
[10, 2, -10, 6],
[6, -14, 10, -12],
[9, 17, -19, 13],
[12, -2, -4, -4],
[5, -11, -4, -9],
[4, 4, 2, 10],
[-27, -3, 5, -7],
[-4, 2, 11, 8],
[-13, -1, -13, 3],
[8, 2, 6, 0],
[-20, 4, 0, -4],
[18, -2, -2, -4],
[-12, 2, 16, -2],
[3, 9, 0, -3],
[21, 3, 3, 13],
[-22, -3, 25, -16],
[-8, 8, -6, 24],
[-17, -15, 13, -1],
[6, 14, 6, 16],
[-2, -8, -13, 10],
[20, -12, -6, -4],
[5, -15, 12, -13]]
