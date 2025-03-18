#for u in GF(q), we can factor as u=aa^*=aa^q=a^{q+1} in GF(q**2) using gen. z and modular arithmetic
def conj_sqrt(u):
    q = u.parent().order().sqrt()
    if u == 0:
        return 0
    z = u.parent().multiplicative_generator()
    k = u.log(z)  # Compute discrete log of u to the base z
    if k % (q+1) != 0:
        raise ValueError("exponent must be divisible by q+1")
    return z ** (k//(q+1))

#given a multiplicative generator `z` of the finite field, the discrete_log is the exponent of the generator
#the discrete_log of zero is -infinity, which we set to -1 for convenience since all other values are nonnegative
def discrete_log(F,x):
    return x.log(F.multiplicative_generator()) if x != 0 else -1

#the map from modular representation theory to compute Brauer characters from F_q --> \C giving roots of unity
#note: Brauer character is the rep'n matrix eigenvalues (over a splitting field extension of F_q) mapped to \C using Brauer map and summed
#i.e. let \alpha = g^k |--> \exp(2*pi*i*k/(q-1))
#i.e. F_q^* is cyclic of order q-1, mapping to (q-1)^th roots of unity in \C
#if the discrete log is provided, we use it directly; otherwise, we compute it
from math import pi
from cmath import exp
def brauer_map(F, a=None, log_a=None):
    """
    Map from F_q to C using the Brauer character formula.
    If `log_a` is provided, it uses that directly; otherwise, it computes the log.
    
    a: Element of F_q or None if log_a is given
    log_a: Precomputed log value (optional)
    F: Finite field F_q
    
    Returns: complex valued root of unity
    """
    if a is None and log_a is None:
        raise ValueError("Either 'a' or 'log_a' must be provided.")
    
    if a is not None:
        if a == 0:
            return 0
        log_a = a.log(F.multiplicative_generator())
    
    return exp(2 * pi * 1j * log_a / (F.order() - 1))

#round a complex number to a given number of digits by rounding real and imaginary parts separately
def round_complex(z, digits):
    try:
        return z.real_part().n(digits=digits) + z.imag_part().n(digits=digits) * 1j
    except AttributeError:
        return z.n(digits=digits)