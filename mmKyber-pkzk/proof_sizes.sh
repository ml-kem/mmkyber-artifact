#!/bin/bash

rm -f proof.log
for level in 128 192 256; do
	make clean
	make LEVEL=$level
	./mm_proof | tee -a proof.log
done
