# discrete-fourier-transform-over-finite-field 
 
Implements the discrete Fourier transform over a ring $R$ in the case where $R$ is a finite field.

The DFT is usually a map $\mathbb{C}^N \rightarrow \mathbb{C}^N$. However, it can be viewed as a decomposition of the group algebra $\mathbb{C}[C_N]$ where $C_N$ is a cyclic group of order $N$. In this case the representations of the cyclic group are one-dimensional $\chi_k(n) = exp(2\pi i k / n)$, the characters. They span the space of functions on $C_N$, and are orthogonal and linearly independent. The Fourier transform is just $x \mapsto \sum_n x_n \chi_k(n)$. This is general for finite groups, where instead of characters we use representations when the group is nonabelian.

We can look at an analogous case over a finite field, the decomposition of $F_p[C_N] = F_p[x]/(x^N-1)$. We have two cases:

$p \nmid N$: The polynomial $x^N-1$ is seperable, so is a product of irreducible factors of multiplicity 1. When $N|p-1$, there exists a primitive root of unity $\alpha$, so we can just use the formula $\displaystyle x \mapsto \sum_{n=0}^{N-1} x_n \alpha^{nk}$ where $\alpha$ is a primitive $N^{th}$ root of unity. This corresponds to the factorization $x^N-1 = (x - \alpha^0) \ldots (x - \alpha^{N-1})$. When $N \nmid p-1$, we have the decomposition $F_p[x]/(x^N-1) = \bigoplus_{d|m, i} F_p[x]/(P_i(x))$ where we factor $x^N-1$ into cyclotomic polynomials $\Phi_d(x)$. We further factor those into $P_1 \ldots P_g$, where each factor has residue degree $f$, and there are $g$ polynomials. Here $fg = \phi(d)$, where $\phi$ is Euler's totient function, and $f$ is the order of $p$ modulo $d$. This is equivalent to factoring the prime ideal (p) in $\mathbb{Z}[\zeta] = \mathbb{Z}[x]/(\Phi_d(x))$

$p | N$: Here we need to consider that $1/N$ does not exist. We factor $x^N-1 = (x^m-1)^{p^s}$. We now repeat the procedure above with $x^m-1$, obtaining a product of cyclotomic polynomials, and then further factor those using algebraic number theory to obtain $F_p[x]/(x^N-1) = \bigoplus_{d|m, i} F_p[x]/(P_i(x)^{p^s})$.

We can then use the Chinese remainder theorem to obtain a linear isomorphism $F_p^N \rightarrow F_p^N$.

We can always find a splitting field $F_q$ of $x^N-1$ for some $q=p^r$. This is equivalent to finding $N|q-1$ for some $r$. This gives a primitive $N^{th}$ root in the first case, splitting completely. When $p|N$, each linear factor has multiplicity $p^s$.

---

For non-abelian groups, representation theory in higher dimensions is required.

For the symmetric group, note that when $p=char(k)$ does not divide $N$, we may use Maschke's theorem to decompose $k[S_N] = \bigoplus_i End(V_i)$ where $V_i$ are the irreducible representations. For the symmetric group, these are labeled by partitions $\lambda$, and $V_i = S^\lambda$ are Specht modules.

When $p|N$, the situation is more complicated and requires modular representation theory. Here the Specht modules are no longer irreducible and Maschke's theorem does not apply. Instead, the simple modules $D^\lambda = S^\lambda / S^\lambda \cap (S^\lambda)^\perp$ are indexed by p-regular partitions of $\lambda$ which do not have any parts with multiplicity more than p, e.g. $(6^2,5,1^3)$ is not 2-regular, since 1 is repeated 3 times. 

We now decompose the group algebra into blocks, $k[S_N] = \bigoplus_i k[S_N]e_i$, where $e_i$ are primitive central orthogonal idempotents. These may be constructed as in Murphy, "The ldempotents of the Symmetric Group and Nakayama’s Conjecture". The modules $D^\lambda$ occur as simple factors in the composition series. The decomposition matrices $D^{\lambda\mu}$ record the number of times the simple module $D^\mu$ occurs in the composition series for $S^\lambda$. The blocks are labeled by p-cores $\gamma$ which are partitions where all possible rim $p$-hooks have been removed. 
