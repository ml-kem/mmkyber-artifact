#!/bin/bash

echo -n "# " > bench.txt
uname -a >> bench.txt
echo -n "# " >> bench.txt
date >> bench.txt

for level in 128 192 256; do
	make obj-clean
	make LEVEL=$level
	echo >> bench.txt
	./xtest | tee -a bench.txt
done

