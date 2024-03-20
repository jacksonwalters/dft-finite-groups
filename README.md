# discrete-fourier-transform-over-finite-field 
 
Implements the discrete Fourier transform over a ring $R$ in the case where $R$ is a finite field.

The DFT is usually a map $\mathbb{C}^N \rightarrow \mathbb{C}^N$. However, it can be viewed as a decomposition of the group algebra $\mathbb{C}[C_N]$ where $C_N$ is a cyclic group of order $N$. In this case the representations of the cyclic group are one-dimensional $\chi_k(n) = exp(2\pi i k / n)$, the characters. They span the space of functions on $C_N$, and are orthogonal and linearly independent. The Fourier transform is just $x \mapsto \sum_n x_n \chi(n)$. This is general for finite groups, where instead of characters we use representations when the group is nonabelian.

We can look at an analogous case over a finite field, the decomposition of $F_p[C_N] = F_p[x]/(x^N-1)$. We have two cases:

$p \nmid N$: When $p$ does not divide $N$, we are just looking at the cyclotomic field of order $N$. All the representations of $C_N$ over $F_p$ are one dimensional and irreducible, and we can just use the formula $x \mapsto \sum_{n=0}^{N-1} x_n \alpha^{nk}$ where $\alpha$ is a primitive $Nth root of unity$.

$p | N$: Here we need to consider that $1/N$ does not exist. Still, we have the decomposition $F_p[x]/(x^N-1) = \prod_{d|m, i} F_p[x]/(P_i(x)^{p^s})$ where we factor (x^N-1)=(x^m-1)^{p^s}. We then factor x^m-1 into cyclotomic polynomials phi_d(x). We further factor those into P_1 ... P_g, where each factor has residue degree f, and there are g polynomials. This is equivalent to factoring the prime ideal (p) in Z[\zeta] = Z[x]/(phi_d(x)). 

Though it is possible to use the Chinese remainder theorem to obtain a linear isomorphism $F_p^N \rightarrow F_p^N$, it may be better to find out how to use a primitive m-th root of unity, and deal with the p^s multiplicity to align this case with the first case.

---

For non-abelian groups, the situation is more complicated and requires modular representation theory.
