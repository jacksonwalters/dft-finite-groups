"""
Constructs the DFT of the dihedral group in SageMath.
Supports both symbolic (complex) and finite field modes.
"""

from sage.all import *

# ── Configuration ──────────────────────────────────────────────────────────────
n = 10; print("n =", n)
USE_FINITE_FIELD = True
q = 11  # prime power; ignored when USE_FINITE_FIELD = False

# ── Build the coefficient ring and choose omega ─────────────────────────────────
if USE_FINITE_FIELD:
    # Find the smallest k such that n | (q^k - 1),
    # i.e. the multiplicative order of q mod n.
    k = Zmod(n)(q).multiplicative_order()
    F = GF(q**k, 'a')
    print(f"Working in GF({q}^{k}) = GF({q**k})")

    # Find a primitive n-th root of unity in F.
    # The multiplicative group of GF(q^k) is cyclic of order q^k - 1,
    # so a generator g satisfies g^((q^k-1)/n) has order n.
    g = F.multiplicative_generator()
    omega = g**((q**k - 1) // n)
    assert omega**n == F.one(), "omega is not an n-th root of unity"
    assert omega**(n-1) != F.one() or n == 1, "omega is not primitive"
    print(f"omega = {omega}  (order {omega.multiplicative_order()})")
else:
    # Symbolic / complex mode (original behaviour)
    var('z')
    forget()
    assume(z**n == 1)
    omega = z
    F = None
    print(f"omega = {omega}  (symbolic)")

# ── Dihedral group setup ────────────────────────────────────────────────────────
G = DihedralGroup(n); print("G =", G)

r = [g for g in G if g.order() == n][0]
s = [g for g in G if g.order() == 2 and g != r**(n//2)][0]
print("r =", r)
print("s =", s)

def express_in_gens(g):
    """Return (0, k) if g = r^k, or (1, k) if g = s·r^k."""
    for k in range(n):
        if g == r**k:
            return (0, k)
    for k in range(n):
        if g == s * r**k:
            return (1, k)

# ── Representation matrices ─────────────────────────────────────────────────────
def mat(entries, ring=None):
    """Build a matrix over `ring` (or infer from entries)."""
    if ring is not None:
        return matrix(ring, entries)
    return matrix(entries)

def rho_odd(k, g):
    (s_exp, r_exp) = express_in_gens(g)
    if k == 0:
        return mat([[F.one() if F else 1]])
    if k == -1:
        return mat([[F.one() if F else 1]]) if s_exp == 0 else mat([[(-F.one()) if F else -1]])
    # k >= 1: 2-dimensional
    wr = omega**(k * r_exp)
    wri = omega**(-k * r_exp)  # = wr^{-1} since omega is a unit
    if s_exp == 0:
        return mat([[wr, 0], [0, wri]])
    else:
        return mat([[0, wr], [wri, 0]])

def rho_even(k, g):
    (s_exp, r_exp) = express_in_gens(g)
    one  =  F.one() if F else 1
    mone = -F.one() if F else -1
    if k == 0:   return mat([[one]])                            # trivial
    if k == -1:  return mat([[one if r_exp % 2 == 0 else mone]])  # det of rotation
    if k == -2:  return mat([[one if s_exp == 0    else mone]])   # det of reflection
    if k == -3:  return mat([[one if (r_exp+s_exp) % 2 == 0 else mone]])  # total sign
    # k >= 1: 2-dimensional
    wr  = omega**(k * r_exp)
    wri = omega**(-k * r_exp)
    if s_exp == 0:
        return mat([[wr, 0], [0, wri]])
    else:
        return mat([[0, wr], [wri, 0]])

# ── Assemble DFT matrix ─────────────────────────────────────────────────────────
def dft_matrix_odd():
    assert n % 2 == 1
    rows = []
    for g in G:
        row = [rho_odd(k, g).list() for k in range(-1, (n-1)//2 + 1)]
        rows.append(sum(row, []))
    return matrix(F, rows) if F else matrix(rows)

def dft_matrix_even():
    assert n % 2 == 0
    rows = []
    for g in G:
        row = [rho_even(k, g).list() for k in range(-3, n//2)]
        rows.append(sum(row, []))
    return matrix(F, rows) if F else matrix(rows)

DFT = dft_matrix_odd() if n % 2 == 1 else dft_matrix_even()

print("\nDFT matrix:")
print(DFT)
print("\nTrace of DFT:", DFT.trace())

# ── Sanity check: columns should be orthogonal (scaled) ────────────────────────
if USE_FINITE_FIELD:
    check = DFT * DFT.transpose()
    is_scalar = all(
        check[i, j] == (F.zero() if i != j else check[i, i])
        for i in range(check.nrows()) for j in range(check.ncols())
    )
    print("\nDFT · DFTᵀ is scalar-block diagonal:", is_scalar)