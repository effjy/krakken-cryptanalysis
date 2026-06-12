# Makefile for Krakken-2048 diffusion tests (SAC / completeness / monobit)
#
#   make            build the sac_test binary
#   make run        build and run the default probe (rounds=8, N=4096)
#   make scan       build and run all rounds 1..8
#   make quick      fast smoke test (rounds=8, N=512)
#   make clean      remove build artifacts
#
# The AVX2 permutation requires a CPU with AVX2 (-mavx2).

CC      ?= gcc
CFLAGS  ?= -O3 -mavx2 -Wall -Wextra
LDLIBS  := -lpthread -lm

BIN     := sac_test
SRCS    := sac_test.c krakken_multi.c

.PHONY: all run scan quick clean

all: $(BIN)

$(BIN): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $@ $(LDLIBS)

run: $(BIN)
	./$(BIN)

scan: $(BIN)
	./$(BIN) all

quick: $(BIN)
	./$(BIN) 8 512

clean:
	rm -f $(BIN)
