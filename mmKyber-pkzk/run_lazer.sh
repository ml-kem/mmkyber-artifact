#!/bin/bash

THIS_DIR=`pwd`
LAZER_DIR="../../lazer"

for pyf in mm_*param.py; do
	hdf=`basename -s .py $pyf`.h
	echo "=== $pyf -> $hdf"
	(cd $LAZER_DIR/scripts && sage lin-codegen.sage \
		$THIS_DIR/$pyf > $THIS_DIR/$hdf)
done
