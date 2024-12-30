#!/bin/bash
cat ../bench-ref.txt ../../mmKyber-c/bench-mm.txt | python3 make_plots.py bytes
cat ../bench-ref.txt ../../mmKyber-c/bench-mm.txt | python3 make_plots.py speed

