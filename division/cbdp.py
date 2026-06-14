"""Conventional bit-based division property (CBDP) as a SCIP MILP.
Primitives: COPY, XOR, linear layer (from GF(2) matrix), S-box (trail selectors).
Rules (Xiang et al. 2016, two-subset CBDP):
  COPY  a -> (b_1..b_m):  sum_i b_i = a
  XOR   (a_1..a_m) -> b:  b = sum_i a_i
  SBOX  u -> v: pick exactly one valid (minimal) division trail.
"""
import pyscipopt as scip
from sbox_trails import all_trails

class CBDP:
    def __init__(self):
        self.m = scip.Model()
        self.m.setParam("display/verblevel", 0)
        self._n = 0
    def bit(self):
        v = self.m.addVar(vtype="B", name=f"x{self._n}")
        self._n += 1
        return v
    def bits(self, k):
        return [self.bit() for _ in range(k)]
    def copy(self, a, m):
        """fan a into m branches: returns m new bits with sum = a."""
        outs = self.bits(m)
        self.m.addCons(scip.quicksum(outs) == a)
        return outs
    def xor(self, ins):
        """returns new bit = division-XOR of ins (sum)."""
        b = self.bit()
        self.m.addCons(b == scip.quicksum(ins))
        return b
    def and_(self, ins):
        """division AND: b = OR of ins.  b >= a_i ;  b <= sum a_i."""
        b = self.bit()
        for a in ins:
            self.m.addCons(b >= a)
        self.m.addCons(b <= scip.quicksum(ins))
        return b
    def const0(self):
        z = self.bit(); self.m.addCons(z == 0); return z
    def linear_layer(self, invars, cols):
        """invars: list of bit vars (len n_in). cols[j]=set of output rows for input j.
        Returns outvars (len n_out). Exact CBDP model: COPY each input over its column
        support, then XOR the routed copies per output row."""
        n_out = 1 + max((max(c) for c in cols if c), default=-1)
        # accumulate, per output row, the list of routed-copy bits
        contrib = {i: [] for i in range(n_out)}
        for j, a in enumerate(invars):
            rows = sorted(cols[j])
            if not rows:
                # input feeds nothing -> must be 0 (can't carry a division bit nowhere)
                self.m.addCons(a == 0)
                continue
            branches = self.copy(a, len(rows))
            for r, c in zip(rows, branches):
                contrib[r].append(c)
        outs = []
        for i in range(n_out):
            if contrib[i]:
                outs.append(self.xor(contrib[i]))
            else:
                z = self.bit(); self.m.addCons(z == 0); outs.append(z)
        return outs
    def sbox(self, invars, trails, n_out=8):
        """invars: 8 input bits. trails: dict u->list of minimal v. Selector encoding."""
        # collect all valid (u,v) pairs
        pairs = [(u, v) for u, vs in trails.items() for v in vs]
        sel = {t: self.bit() for t in range(len(pairs))}
        self.m.addCons(scip.quicksum(sel.values()) == 1)   # exactly one trail
        outs = self.bits(n_out)
        for bit in range(8):
            self.m.addCons(invars[bit] ==
                scip.quicksum(sel[t] for t,(u,v) in enumerate(pairs) if (u>>bit)&1))
        for bit in range(n_out):
            self.m.addCons(outs[bit] ==
                scip.quicksum(sel[t] for t,(u,v) in enumerate(pairs) if (v>>bit)&1))
        return outs
    def sbox_ineq(self, invars, facets, n_out=8):
        """compact facet S-box model: NO extra vars, just inequalities over the
        16 bit-vars [in0..in7, out0..out7].  facets: list of (a[16], b) meaning a.x >= b."""
        outs = self.bits(n_out)
        x = list(invars) + outs
        for a, b in facets:
            self.m.addCons(scip.quicksum(a[i]*x[i] for i in range(16) if a[i]) >= b)
        return outs
    def fix(self, var, val):
        self.m.addCons(var == val)
    def feasible(self):
        self.m.optimize()
        return self.m.getStatus() == "optimal"

# ---------------- correctness gate: MILP S-box must match the trail table ----------------
if __name__ == "__main__":
    T = all_trails()
    # Reconstruct, via the MILP, the set of v reachable from each u; compare to T.
    # For a handful of u, enumerate reachable v by solving feasibility for each candidate v.
    import itertools
    def popc(x): return bin(x).count("1")
    def milp_reach(u):
        found = set()
        for v in range(256):
            c = CBDP()
            ins = c.bits(8); outs = c.sbox(ins, T)
            for b in range(8):
                c.fix(ins[b], (u>>b)&1)
                c.fix(outs[b], (v>>b)&1)
            if c.feasible():
                found.add(v)
        # reduce to minimal
        return sorted(v for v in found if not any(v2!=v and (v2&~v)==0 for v2 in found))
    ok = True
    for u in [0x00, 0x01, 0x03, 0x0F, 0x7F, 0xFF, 0x55, 0x80]:
        got = milp_reach(u); exp = T[u]
        match = (got == exp)
        ok &= match
        print(f"u={u:#04x}: MILP minimal-reach == trail table? {match}  (|exp|={len(exp)})")
    # degree check via MILP: max weight u reaching a single output bit
    degs=[]
    for bit in range(8):
        mx=0
        for u in range(256):
            c=CBDP(); ins=c.bits(8); outs=c.sbox(ins,T)
            for b in range(8): c.fix(ins[b],(u>>b)&1)
            for b in range(8): c.fix(outs[b], 1 if b==bit else 0)
            if c.feasible(): mx=max(mx,popc(u))
        degs.append(mx)
    print("per-bit degree via MILP:", degs, "(expected all 7)")
    print("GATE PASSED" if ok and degs==[7]*8 else "GATE FAILED")
