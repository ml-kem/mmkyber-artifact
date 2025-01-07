# Lattice Parameter Selection for mmKyber

These SageMath scripts implement the exploration of lattice parameters for mmKyber.

The scripts allow users to quickly select the basic parameters, including the modulus $q$, the degree $d$, the dimensions $m$, $n$, the infinity-norm bound $\nu$ and its support $\bar{\nu}$, the Gaussian widths $\sigma_0$, $\sigma_1$, and the compression parameters $d_v$, $d_u$, for different recipient numbers $N$ and security levels of $128$-bit, $192$-bit, or $256$-bit.

Some parameters in the reference implementation have been further tuned manually:

* The compression parameter $d_u$ has been further reduced by performing the decryption failure test in `pr-fail-dec`, leading to slightly shorter ciphertexts.

* The modulus $q=2^{25}-2^{12}+1 = 33550337$ used by the implementation allows more options for NTT than $2^{25}-5 \cdot 2^{8}+1 = 33553153$ output by the script.

