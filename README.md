# DFTs of finite groups over various fields
 
Implements various DFTs for the cyclic group and symmetric group. We work over the complex numbers, finite fields, and number fields.

For finite fields, when $p$ divides $|G|$ we use the modular DFT of the symmetric group, which is given by the Peirce decomposition using central primitive orthogonal idempotents.

For finite fields, when $p > n$ we can form the unitary DFT of the symmetric group. We compute the usual seminormal DFT, multiply by its conjugate transpose, and then factor the resulting diagonal matrix using conjugate square roots, and multiply by the inverse. This normalizes the matrix to be unitary with respect to finite field conjugation.

For unitary representations of the symmetric group over finite fields, we compute an $S_n$-invariant symmetric bilinear form using linear algebra, then use the extended Cholesky decomposition to factor the matrix associated to the form. This yields a change-of-basis matrix which makes the representation unitary.

Over number fields, we follow the same pattern and multiply by the conjugate transpose to get a diagonal matrix. Then we compute a number field which contains all square roots of those elements, and then take square roots along the diagonal, and multiply by the inverse of the resulting diagonal matrix.

See:  [The Modular DFT of the Symmetric Group](https://doi.org/10.48550/arXiv.2404.05796)
