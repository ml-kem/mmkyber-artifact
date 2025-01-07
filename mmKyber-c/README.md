#   mmKyber-c

Proof-of-Concept implementation of mmKyber-KEM and mmKyber-PKE in plain C.


##  C Implementation Notes

Plain C code, no library dependencies or assembly optimization.
However -- currently only tested with an x86 Linux.

Some components (especially Gaussian samplers) are temporary.
Furthermore this code is not consistently constant-time, although
attempt has been made in some places.


##  Reproducing benchmarks

A file `bench.txt` is produced by `make bench`.


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

