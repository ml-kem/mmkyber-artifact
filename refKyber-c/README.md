#	refKyber-c

This is the benchmark code used for "apples-to-apples" comparison between
plain C implementations of Kyber and mmKyber.

There is a similar `Makefile` to the mmKyber benchmark, with equivalent
compiler settings. The script names are also the same; `run_bench.sh`
creates a text file `bench.txt`.
```
./run_bench.sh
```

For purposes of reporting, the function names in the Kyber implementation
have been mapped to the ML-KEM (CCA) and K-PKE (CPA) names from 
[FIPS 203](https://doi.org/10.6028/NIST.FIPS.203).


##	Reference Kyber code

The `kyber` symbolic link should point to the official 
[reference implementation](https://github.com/pq-crystals/kyber) of
Kyber from the Kyber design team. This code was tested against the
commit `10b478f` from August 21, 2024 (after FIPS 203 came out.)

You can also have the repo right here:
```
rm -f kyber
git clone https://github.com/pq-crystals/kyber.git
```

