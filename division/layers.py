"""Probe Krakken layers via ctypes; extract GF(2) matrices for linear layers
and verify linearity. State = 32 uint64 words = 2048 bits.
Bit index convention: global bit b -> word b//64, bit-in-word b%64 (LSB-first).
"""
import ctypes, os, random
NB = 2048
LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "libkrakken.so")
_lib = ctypes.CDLL(LIB)
StateT = ctypes.c_uint64 * 32

LINEAR_LAYERS = ["theta_scalar","tentacle_mds_scalar","rho_scalar","pi_scalar",
                 "butterfly_diffusion_scalar","ink_cloud_shuffle"]
NONLINEAR     = ["chi_scalar","pressure_arx_scalar"]  # sbox / ARX

def _apply(fn, bits):
    """bits: iterable of set bit indices -> returns set of output bit indices."""
    st = StateT()
    for b in bits:
        st[b >> 6] |= (1 << (b & 63))
    fn(st)
    out = set()
    for w in range(32):
        v = st[w]
        while v:
            lsb = v & (-v)
            out.add(w*64 + lsb.bit_length()-1)
            v ^= lsb
    return out

def extract_matrix(name):
    """Return columns: matrix[j] = set of output bits set by input unit vector e_j.
    (i.e. column j of the GF(2) matrix). Assumes the layer is linear with F(0)=0."""
    fn = getattr(_lib, name)
    cols = []
    for j in range(NB):
        cols.append(frozenset(_apply(fn, [j])))
    return cols

def verify_linear(name, cols, trials=40):
    fn = getattr(_lib, name)
    # F(0) == 0 ?
    if _apply(fn, []):
        return False, "F(0) != 0"
    rng = random.Random(12345)
    for _ in range(trials):
        k = rng.randint(1, 8)
        A = set(rng.sample(range(NB), k))
        B = set(rng.sample(range(NB), k))
        fa = _apply(fn, A); fb = _apply(fn, B)
        fab = _apply(fn, A ^ B)            # symmetric difference = XOR of inputs
        # predicted from columns: XOR of columns over set bits
        pred = set()
        for j in (A ^ B):
            pred ^= set(cols[j])
        if fab != (fa ^ fb):
            return False, "not additive on random sample"
        if fab != pred:
            return False, "column-sum mismatch"
    return True, "ok"

if __name__ == "__main__":
    for name in LINEAR_LAYERS:
        cols = extract_matrix(name)
        ok, msg = verify_linear(name, cols)
        nz = sum(len(c) for c in cols)
        # invertibility quick check: are all columns distinct & none empty?
        empties = sum(1 for c in cols if not c)
        print(f"{name:28s} linear={ok!s:5s} ({msg}); nnz/col avg={nz/NB:5.2f}; empty cols={empties}")
    for name in NONLINEAR:
        fn = getattr(_lib, name)
        zero = _apply(fn, [])  # F(0)
        # additivity test should FAIL for these
        rng = random.Random(7)
        A={1,2}; B={3,4}
        add = _apply(fn,A^B)==(_apply(fn,A)^_apply(fn,B))
        print(f"{name:28s} NONLINEAR (F(0) set bits={len(zero)}, additive_on_sample={add})")
