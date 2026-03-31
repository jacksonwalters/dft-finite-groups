"""
Constructs the DFT of the dihedral group in SageMath using knowledge of the representation theory.
""";

from sage.all import *

#construct dihedral group of order n over GF(q)
n = 10; print("n =", n)
var('z')
forget()
assume(z**n == 1)
omega = z; print("omega =", omega)
G = DihedralGroup(n); print("G =", G)
#gens = G.gens()
r = [g for g in G if g.order() == n][0]
s = [g for g in G if g.order() == 2 and g != r**(n//2)][0]
print("r =", r)
print("s =", s)

# returns (0, k) if g = r^k and (1, k) if g = s*r^k
def express_in_gens(g):
    for k in range(n):
        if g == r**k:
            return (0, k)
    for k in range(n):
        if g == s * r**k:
            return (1, k)
        
# n odd, we have two 1-dim'l irreps and (n-1)/2 2-dim'l irreps
# the 1-dim's irreps are trivial and sign
# the 2-dim'l irreps are given by rotation matrices and a flip matrix
def rho_odd(k, g):
    (s_exp, r_exp) = express_in_gens(g)
    if k == 0:
        return matrix([1])
    if k == -1:
        if s_exp == 0:
            return matrix([1])
        if s_exp == 1:
            return matrix([-1])
    if k >= 1:
        if s_exp == 0:
            return matrix([[omega**(k*r_exp), 0], [0, omega**(-k*r_exp)]])
        if s_exp == 1:
            return matrix([[0, omega**(k*r_exp)], [omega**(-k*r_exp), 0]])
        
def dft_matrix_odd():
    assert n % 2 == 1
    rows = []
    for g in G:
        row = [rho_odd(k, g).list() for k in range(-1,(n-1)//2 + 1)]
        rows.append(sum(row, []))
    return matrix(rows)

# for n even case
def rho_even(k, g):
    (s_exp, r_exp) = express_in_gens(g)
    if k == 0:   # trivial
        return matrix([1])
    if k == -1:  # sign of rotation
        return matrix([(-1)**r_exp])
    if k == -2:  # sign of reflection
        return matrix([(-1)**s_exp])
    if k == -3:  # total sign
        return matrix([(-1)**(r_exp + s_exp)])
    if k >= 1:
        if s_exp == 0:
            return matrix([[omega**(k*r_exp), 0], [0, omega**(-k*r_exp)]])
        if s_exp == 1:
            return matrix([[0, omega**(k*r_exp)], [omega**(-k*r_exp), 0]])
        
# form the DFT matrix for n even
def dft_matrix_even():
    assert n % 2 == 0
    rows = []
    for g in G:
        row = [rho_even(k, g).list() for k in range(-3, n//2)]
        rows.append(sum(row, []))
    return matrix(rows)

DFT_matrix = dft_matrix_odd() if n % 2 == 1 else dft_matrix_even(); 

print(DFT_matrix)

print("trace of DFT:", DFT_matrix.trace())