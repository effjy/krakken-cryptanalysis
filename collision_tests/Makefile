# Makefile for Krakken-2048 sponge-hash collision / distribution battery
#
#   make            build the collision_test binary
#   make run        build and run with default N = 1,000,000
#   make quick      fast smoke test (N = 100,000)
#   make clean      remove build artifacts
#
# Requires a CPU with AVX2 (-mavx2).

CC      ?= gcc
CFLAGS  ?= -O3 -mavx2 -Wall -Wextra
LDLIBS  := -lpthread -lm

BIN     := collision_test
SRCS    := collision_test.c krakken_multi.c

.PHONY: all run quick clean

all: $(BIN)

$(BIN): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $@ $(LDLIBS)

run: $(BIN)
	./$(BIN)

quick: $(BIN)
	./$(BIN) 100000

clean:
	rm -f $(BIN)
