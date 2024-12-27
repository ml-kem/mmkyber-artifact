#!/bin/bash
cat ../bench.txt ../../mmKyber-c/bench.txt | python3 make_plots.py bytes
cat ../bench.txt ../../mmKyber-c/bench.txt | python3 make_plots.py speed

