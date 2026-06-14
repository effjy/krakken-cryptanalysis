"""Krakken round model for CBDP. chi & pressure_arx are defined as abstract
gate circuits with two interpreters (concrete / division) so the wiring can be
validated against the C reference and reused unchanged in the MILP."""
import ctypes, os, pickle, functools
from sbox_trails import ABYSSAL, all_trails
import layers as L

HERE = os.path.dirname(os.path.abspath(__file__))
NB = 2048

# ---------- abstract circuit engines ----------
class Concrete:
    """plain-bit interpreter (for wiring validation)."""
    def copy(self, a, m): return [a]*m
    def xor(self, ins):
        r=0
        for x in ins: r ^= x
        return r
    def and_(self, ins):
        r=1
        for x in ins: r &= x
        return r
    def const0(self): return 0
    def sbox(self, in8):
        val=0
        for j,b in enumerate(in8): val |= (b&1)<<j
        o=ABYSSAL[val]
        return [(o>>j)&1 for j in range(8)]

class DivEngine:
    """MILP interpreter wrapping a CBDP model. Uses compact facet S-box model if
    facets are available (facets.pkl), else the per-trail selector encoding."""
    def __init__(self, c, trails, facets=None):
        self.c=c; self.trails=trails; self.facets=facets
    def copy(self, a, m): return self.c.copy(a, m)
    def xor(self, ins): return self.c.xor(list(ins))
    def and_(self, ins): return self.c.and_(list(ins))
    def const0(self): return self.c.const0()
    def sbox(self, in8):
        if self.facets is not None:
            return self.c.sbox_ineq(list(in8), self.facets, n_out=8)
        return self.c.sbox(list(in8), self.trails, n_out=8)

# ---------- word bit-list helpers (index i = bit i, LSB=0) ----------
def rotl(lst, r): return [lst[(p-r)%64] for p in range(64)]
def shr(eng, lst, s): return [lst[p+s] if p+s<64 else eng.const0() for p in range(64)]
def shl(eng, lst, s): return [lst[p-s] if p-s>=0 else eng.const0() for p in range(64)]
def xorw(eng, x, y): return [eng.xor([a,b]) for a,b in zip(x,y)]

def sbox64(eng, w):
    out=[None]*64
    for k in range(8):
        o=eng.sbox(w[8*k:8*k+8])
        for j in range(8): out[8*k+j]=o[j]
    return out

# ---------- chi circuit (one column-pair) ----------
def chi_pair(eng, a, b):
    # b used in t1 (rotl 32) and t2 -> split each bit into 2
    bc=[eng.copy(bit,2) for bit in b]; b1=[p[0] for p in bc]; b2=[p[1] for p in bc]
    t1=xorw(eng, a, rotl(b1,32))      # a used once
    ap=sbox64(eng, t1)
    apc=[eng.copy(bit,2) for bit in ap]; ap1=[p[0] for p in apc]; ap2=[p[1] for p in apc]
    t2=xorw(eng, b2, rotl(ap2,32))
    bp=sbox64(eng, t2)
    return ap1, bp                     # -> word ca, word cb

# ---------- modular add (mod 2^64) exact full-adder division circuit ----------
def add64(eng, x, y):
    out=[]; carry=eng.const0()
    for i in range(64):
        xs=eng.copy(x[i],3); ys=eng.copy(y[i],3); cs=eng.copy(carry,3)
        s=eng.xor([xs[0], ys[0], cs[0]])
        m1=eng.and_([xs[1], ys[1]]); m2=eng.and_([xs[2], cs[1]]); m3=eng.and_([ys[2], cs[2]])
        carry=eng.xor([m1,m2,m3])
        out.append(s)
    return out                         # carry out of bit 63 discarded (mod 2^64)

# ---------- pressure_arx circuit (one column: words a,b,cc,d) ----------
def arx_col(eng, a, b, cc, d):
    cc3=[eng.copy(bit,3) for bit in cc]; ccA=[p[0] for p in cc3]; ccB=[p[1] for p in cc3]; ccBase=[p[2] for p in cc3]
    d3 =[eng.copy(bit,3) for bit in d ]; dA =[p[0] for p in d3 ]; dB =[p[1] for p in d3 ]; dBase=[p[2] for p in d3]
    g_a=xorw(eng, ccA, shr(eng, ccB,17))
    g_b=xorw(eng, dA,  shr(eng, dB, 17))
    a1=add64(eng, a, g_a)
    b1=add64(eng, b, g_b)
    a1c=[eng.copy(bit,3) for bit in a1]; a1o=[p[0] for p in a1c]; a1d=[p[1] for p in a1c]; a1s=[p[2] for p in a1c]
    b1c=[eng.copy(bit,3) for bit in b1]; b1o=[p[0] for p in b1c]; b1d=[p[1] for p in b1c]; b1s=[p[2] for p in b1c]
    cc1=add64(eng, ccBase, xorw(eng, a1d, shl(eng, a1s,31)))
    d1 =add64(eng, dBase,  xorw(eng, b1d, shl(eng, b1s,31)))
    return a1o, rotl(b1o,7), cc1, rotl(d1,19)   # outputs a,b,cc,d

# ---------- linear-layer matrices (probe & cache) ----------
def _probe_seq(funcs):
    """compose a sequence of linear C layers -> cols[gid_in] = set(gid_out)."""
    fns=[getattr(L._lib, n) for n in funcs]
    cols=[]
    for j in range(NB):
        st=L.StateT(); st[j>>6] |= (1<<(j&63))
        for fn in fns: fn(st)
        s=set()
        for w in range(32):
            v=st[w]
            while v:
                lsb=v&-v; s.add(w*64+lsb.bit_length()-1); v^=lsb
        cols.append(frozenset(s))
    return cols

@functools.lru_cache(maxsize=None)
def matrices():
    cache=os.path.join(HERE,"matrices.pkl")
    if os.path.exists(cache):
        with open(cache,"rb") as f: return pickle.load(f)
    M={"L1":_probe_seq(["theta_scalar","tentacle_mds_scalar","rho_scalar","pi_scalar"]),
       "BUT":_probe_seq(["butterfly_diffusion_scalar"]),
       "INK":_probe_seq(["ink_cloud_shuffle"])}
    with open(cache,"wb") as f: pickle.dump(M,f)
    return M

# ---------- full round on a div model ----------
CHI_PAIRS=[(((2*p)*4+y),((2*p+1)*4+y)) for y in range(4) for p in range(4)]
ARX_COLS=[(4*c,4*c+1,4*c+2,4*c+3) for c in range(8)]

def state_words(state):  # state: 2048 vars -> 32 lists of 64
    return [state[w*64:w*64+64] for w in range(32)]
def flatten(words):
    out=[None]*NB
    for w in range(32):
        for i in range(64): out[w*64+i]=words[w][i]
    return out

def chi_layer(eng, state):
    words=state_words(state); newords=[None]*32
    for ca,cb in CHI_PAIRS:
        ap,bp=chi_pair(eng, words[ca], words[cb])
        newords[ca]=ap; newords[cb]=bp
    return flatten(newords)

def arx_layer(eng, state):
    words=state_words(state); newords=[None]*32
    for a,b,cc,d in ARX_COLS:
        wa,wb,wcc,wd=arx_col(eng, words[a],words[b],words[cc],words[d])
        newords[a]=wa; newords[b]=wb; newords[cc]=wcc; newords[d]=wd
    return flatten(newords)

def build_rounds(c, trails, state, rounds, verbose=True, facets=None):
    """apply `rounds` Krakken rounds to division state using CBDP model c."""
    import time
    M=matrices(); eng=DivEngine(c,trails,facets)
    for r in range(rounds):
        t=time.time()
        state=c.linear_layer(state, M["L1"])     # theta,mds,rho,pi
        state=chi_layer(eng, state)              # nonlinear S-box layer
        state=c.linear_layer(state, M["BUT"])    # butterfly
        state=arx_layer(eng, state)              # nonlinear ARX layer
        state=c.linear_layer(state, M["INK"])    # ink_cloud (beta_iota const ignored)
        if verbose:
            print(f"    [build] round {r+1}/{rounds} done "
                  f"({c.m.getNVars()} vars, {c.m.getNConss()} cons, "
                  f"+{time.time()-t:.1f}s)", flush=True)
    return state
