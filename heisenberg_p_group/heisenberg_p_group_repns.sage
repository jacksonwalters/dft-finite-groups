from sage.all import *

p = 3
F = GF(p)

G = [(a,b,c) for a in F for b in F for c in F]

K.<zeta> = CyclotomicField(p)
R.<x> = PolynomialRing(K)
L = K.extension(x^2 - len(G), 'a')

def rho_k(k):
    z = zeta**k
    def rep(a,b,c):
        M = matrix(K, p, p)
        for x in range(p):
            y = (x + int(a)) % p
            M[y, x] = z**(int(c) + int(b)*x)
        return M
    return rep

reps = []
for k in F:
    if k != 0:
        reps.append(rho_k(k))

def chi(u, v):
    def rep(a, b, c):
        return zeta**(int(u)*int(a) + int(v)*int(b))
    return rep

def multiply(g1, g2):
    a1,b1,c1 = g1
    a2,b2,c2 = g2
    return (a1+a2, b1+b2, c1+c2 + a1*b2)

u, v = F(1), F(2)
chi_uv = chi(u,v)

g = (F(1),F(2),F(3))
h = (F(2),F(1),F(4))

chi_uv(*multiply(g,h)) == chi_uv(*g) * chi_uv(*h)

def all_irreps():
    irreps = []

    # 1D irreps
    for u in F:
        for v in F:
            irreps.append(("1d", chi(u,v)))

    # p-dimensional irreps
    for k in F:
        if k != 0:
            irreps.append(("pd", rho_k(k)))

    return irreps

def fourier_matrix():
    irreps = all_irreps()
    N = len(G)   # p^3

    rows = []

    for typ, rep in irreps:
        if typ == "1d":
            row = []
            for g in G:
                row.append(rep(*g))
            rows.append(row)

        else:  # p-dimensional
            d = p
            for i in range(d):
                for j in range(d):
                    row = []
                    for g in G:
                        row.append(rep(*g)[i,j])
                    rows.append(row)

    return matrix(K, rows)

def fourier_matrix_normalized():
    irreps = all_irreps()
    N = len(G)

    rows = []

    for typ, rep in irreps:
        if typ == "1d":
            scale = L(1) / L.gen()   # 1 / sqrt(N)
            row = [scale * L(rep(*g)) for g in G]
            rows.append(row)

        else:
            d = p
            scale = L(1) / L.gen()
            for i in range(d):
                for j in range(d):
                    row = [scale * L(rep(*g)[i,j]) for g in G]
                    rows.append(row)

    return matrix(L, rows)

DFT = fourier_matrix(); 

print("DFT eigenvalues:", DFT.eigenvalues())

uDFT = fourier_matrix_normalized()