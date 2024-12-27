#!/bin/bash

rm -f testvec.tmp
for level in 128 192 256; do
	for mode in KEM PKE; do
		make obj-clean
		make MODE=$mode LEVEL=$level TEST=TESTVEC
		./xtest | grep chk >> testvec.tmp
	done
done
sha256sum testvec.tmp

