#	Makefile

test:
	python3 test_vec.py | grep chk > testvec.tmp
	sha256sum testvec.tmp
	cmp testvec.tmp testvec.txt

clean:
	rm -rf __pycache__ testvec.tmp
	rm	-f *-data-*.h

