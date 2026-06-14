# Krakken-2048 division-property integral search.
#
#   make            build libkrakken.so + generate the exact facet S-box model
#   make validate   run the correctness gates (S-box, wiring, ARX rule)
#   make clean      remove generated artifacts
#
# The Python venv with pyscipopt/SCIP, numpy and scipy is required.
# Override its interpreter with:  make PY=/path/to/python

PY     ?= $(HOME)/venv/bin/python
CC     ?= gcc
CFLAGS ?= -shared -fPIC -O2

.PHONY: all validate clean

all: libkrakken.so facets.pkl

# Reference permutation, probed by layers.py to extract the linear matrices.
libkrakken.so: krakken.c krakken.h
	$(CC) $(CFLAGS) krakken.c -o $@ -lpthread

# Compact, EXACT inequality model for the ABYSSAL S-box division polytope.
# This is what makes the MILP solvable (it removes ~516k selector binaries/round).
facets.pkl: gen_facets.py sbox_trails.py
	$(PY) gen_facets.py --reduce

validate: libkrakken.so
	$(PY) validate_wiring.py
	$(PY) validate_arx_rule.py

clean:
	rm -f libkrakken.so facets.pkl matrices.pkl
	rm -rf __pycache__
