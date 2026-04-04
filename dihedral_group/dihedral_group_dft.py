"""
Constructs the DFT of the dihedral group in SageMath using knowledge of the representation theory.

TO_DO: determine why DFT matrix is not unitary

""";

# ── Configuration ──────────────────────────────────────────────────────────────
n = 5; print("n =", n)
USE_FINITE_FIELD = False
p = 23; print("p =", p) # characteristic of the finite field
q = 23; print("q =", q) # size of the finite field

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
    K.<z> = CyclotomicField(n) #cyclotomic field containing a primitive n-th root of unity
    omega = z; print(f"omega = {omega}  (primitive n-th root of unity)")

G = DihedralGroup(n); print("G =", G) #D_n, dihedral group of order 2n
r = [g for g in G if g.order() == n][0] #rotation of order n, generator
s = [g for g in G if g.order() == 2 and g != r**(n//2)][0] #flip of order 2, generator
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
def rho_odd(k, g, omega):
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
            return matrix([[0, omega**(-k*r_exp)], [omega**(k*r_exp), 0]])
        
def dft_matrix_odd(omega,unitary=False):
    assert n % 2 == 1
    dim = lambda k: 1 if k in (-1, 0) else 2
    rows = []
    for g in G:
        if unitary:
            row = [(sqrt(dim(k)/(2*n))*rho_odd(k, g, omega)).list() for k in range(-1,(n-1)//2 + 1)]
        else:
            row = [rho_odd(k, g, omega).transpose().list() for k in range(-1,(n-1)//2 + 1)]
        rows.append(sum(row, []))
    return matrix(rows)

# for n even case
def rho_even(k, g, omega):
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
            return matrix([[0, omega**(-k*r_exp)], [omega**(k*r_exp), 0]])
        
# form the DFT matrix for n even
def dft_matrix_even(omega):
    assert n % 2 == 0
    rows = []
    for g in G:
        row = [rho_even(k, g, omega).list() for k in range(-3, n//2)]
        rows.append(sum(row, []))
    return matrix(rows)

# compute the trace of the DFT matrix using the formula for the entries
def matrix_entries(i, j):
    if j == 0:
        return 1 if i < n else -1
    if j == 1:
        return 1

    k = j - 2
    m = k // 4 + 1
    r = k % 4

    if i < n:
        return omega**(m * i) if r == 0 else omega**(-m * i) if r == 3 else 0
    else:
        return omega**(m * i) if r == 1 else omega**(-m * i) if r == 2 else 0
    
# compute the trace of the DFT matrix using the formula for the entries
def trace_dft():
    total = 2
    total += sum(omega**(((i-2)//4+1)*i) for i in range(2, n) if i % 4 == 2)
    total += sum(omega**(-((i-2)//4+1)*i) for i in range(2, n) if i % 4 == 1)
    total += sum(omega**(((i-2)//4+1)*i) for i in range(n, 2*n) if i % 4 == 3)
    total += sum(omega**(-(((i-2)//4+1)*i)) for i in range(n, 2*n) if i % 4 == 0)
    return total

def DFT_matrix_odd_from_entries():
    assert n % 2 == 1
    return matrix(2*n, 2*n, matrix_entries)

uDFT = dft_matrix_odd(omega,unitary=True); print(uDFT*uDFT.conjugate_transpose())

DFT_matrix = dft_matrix_odd(omega) if n % 2 == 1 else dft_matrix_even(omega); print(DFT_matrix)

# the determinant is 2*n**n
f = DFT_matrix.charpoly(); print(f)

if USE_FINITE_FIELD:
    L.<a> = f.splitting_field(); print(L)

if USE_FINITE_FIELD:
    R_L.<x> = PolynomialRing(L)
    f_L = R_L(f)
    f_L.factor()

if USE_FINITE_FIELD:
    eigenvalues = f.roots(L, multiplicities=False)
    print(eigenvalues)
    print(len(eigenvalues))

def frobenius(x, p):
    return x**p

def frobenius_orbit(alpha, p):
    orbit = []
    seen = set()
    
    x = alpha
    while x not in seen:
        seen.add(x)
        orbit.append(x)
        x = x^p
    
    return orbit

def frobenius_orbits(eigenvalues, p):
    orbits = []
    seen = set()
    
    for lam in eigenvalues:
        if lam not in seen:
            orb = frobenius_orbit(lam, p)
            orbits.append(orb)
            seen.update(orb)
    
    return orbits

def discrete_log_Fq(x, alpha=None):
    F = x.parent()
    if x == 0:
        raise ValueError("Log undefined for 0")
    if alpha is None:
        alpha = F.multiplicative_generator()
    return x.log(alpha)

if USE_FINITE_FIELD:
    frobenius_orbits(eigenvalues, p)

if USE_FINITE_FIELD:
    len(frobenius_orbit(eigenvalues[0], p))

if USE_FINITE_FIELD:
    dlogs = [discrete_log_Fq(x, alpha=None) for x in eigenvalues]; dlogs

def brauer_map(x):
    if x == 0:
        return 0
    else:
        l = discrete_log_Fq(x, alpha=None)
        return exp(2 * pi * I * l / (L.order()-1))
    
if USE_FINITE_FIELD:
    complex_eigs = [brauer_map(eig) for eig in eigenvalues]; complex_eigs
    eigs_num = [complex(e.n(digits=30)) for e in complex_eigs]; eigs_num

if USE_FINITE_FIELD:
    P = list_plot(
        [(z.real, z.imag) for z in eigs_num],
        plotjoined=False,
        marker='o'
    )

    P.show()

# can normalize by 1/sqrt(D) to get a unitary matrix
if not USE_FINITE_FIELD:
    print(DFT_matrix.conjugate_transpose() * DFT_matrix)