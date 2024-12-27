#   mmKyber-c

PoC implementation of mmKyber-KEM and mmKyber-PKE in plain C.


##  C mplementation characteristics

- Plain C code, no library dependencies or assembly optimization.
However currently tested only with an x86 Linux. The biggest consumer of
cycles is the Keccak f1600 function.

= Some components (especially Gaussian samplers) are temporary;
*work in progress*. Furthermore this code is not consistently constant-time.
The non-constant time decryption/decapsulation function is about 10%
faster than the NTT version on some targets; both are provided
(complile-time flag `MM_NON_CT`.)

Here are benchmarks (cycles and wall-clock seconds) for various parameter
sets with my laptop, AMD Ryzen 7 7840U, compiled with gcc 14.2.0-8 on
Linux 6.12.5-1.

**( Insert table and comparisons.)**


##  Reproducing benchmarks

A file `bench.txt` is produced by `make test`.


##  Functional testing

The test program created by the Makefile is named `xtest`, please see options
in the Makefile. One can pass them through the command line is:

```
make MODE=PKE LEVEL=192 TEST=TESTVEC
./xtest
```

You will have to `make clean` when changing from one parameter set to another,
since the parameters are defined as macros.

There are some known-good checksums in `testvec.txt`, which are in a format
also produced by the Python implementation. Invoking `make test` produces
`testvec.tmp` with all parameter sets, and compares this to `testvec.txt`.

