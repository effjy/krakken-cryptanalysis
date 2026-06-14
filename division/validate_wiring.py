"""Validate chi & pressure_arx CIRCUIT WIRING by running the concrete interpreter
against the real C layers, bit-for-bit on random states."""
import random, ctypes
import model, layers as L
from model import Concrete, chi_layer, arx_layer

def state_to_bits(st):
    return [ (st[w]>>i)&1 for w in range(32) for i in range(64) ]
def bits_to_state(bits):
    st=L.StateT()
    for g,b in enumerate(bits):
        if b: st[g>>6] |= (1<<(g&63))
    return st
def apply_C(fn, bits):
    st=bits_to_state(bits); fn(st); return state_to_bits(st)

eng=Concrete()
rng=random.Random(2024)
def check(name, c_fn, circuit_layer, trials=200):
    bad=0
    for _ in range(trials):
        bits=[rng.randint(0,1) for _ in range(2048)]
        got=circuit_layer(eng, bits)
        exp=apply_C(c_fn, bits)
        if got!=exp: bad+=1
    print(f"{name}: {trials-bad}/{trials} states match C reference" + (" -- OK" if bad==0 else " -- MISMATCH"))
    return bad==0

ok1=check("chi    ", L._lib.chi_scalar, chi_layer)
ok2=check("arx    ", L._lib.pressure_arx_scalar, arx_layer)
print("WIRING OK" if ok1 and ok2 else "WIRING FAILED")
