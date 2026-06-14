#!/usr/bin/env python3
"""Krakken-2048 bit-based division-property integral search.

Finds output bits PROVABLY balanced (zero-sum) over a chosen input cube after a
chosen number of rounds. A bit is balanced iff its unit vector is NOT reachable
by any division trail (MILP infeasible). Sound: any balanced bit reported is a
real integral-distinguisher bit. The model is conservative (esp. ARX), so it may
miss some balanced bits but never reports a false one.

NO internal timeout is imposed. Run as long as you like, or wrap with shell
`timeout`, or pass --solver-time-limit S for a per-bit SCIP cap (a bit hitting
the cap is reported as 'unknown', never as balanced).

Examples:
  python3 krakken_divprop.py --rounds 1 --cube word:0 --max-bits 128
  python3 krakken_divprop.py --rounds 5 --cube dim:30:0 --verify-empirical
"""
import argparse, sys, ctypes, time
from cbdp import CBDP
from sbox_trails import all_trails
import model, layers as L

_T0 = time.time()
def log(msg):
    print(f"[{time.time()-_T0:7.1f}s] {msg}", flush=True)

def parse_cube(spec):
    """return sorted list of active (cube) global bit indices."""
    kind, *rest = spec.split(":")
    if kind == "word":
        w = int(rest[0]); return list(range(w*64, w*64+64))
    if kind == "byte":
        w = int(rest[0]); b = int(rest[1]); base = w*64 + b*8; return list(range(base, base+8))
    if kind == "bits":
        return sorted(int(x) for x in rest[0].split(","))
    if kind == "dim":
        import random
        d = int(rest[0]); seed = int(rest[1]) if len(rest) > 1 else 0
        return sorted(random.Random(seed).sample(range(model.NB), d))
    raise SystemExit(f"bad --cube spec: {spec}")

def build(rounds, cube, trails, facets=None):
    c = CBDP()
    state = c.bits(model.NB)
    active = set(cube)
    for g in range(model.NB):
        c.fix(state[g], 1 if g in active else 0)   # initial division = cube indicator
    out = model.build_rounds(c, trails, state, rounds, facets=facets)
    return c, out

def run(args):
    cube = parse_cube(args.cube)
    log(f"=== Krakken-2048 division-property integral search ===")
    log(f"rounds={args.rounds}  cube='{args.cube}' (dim {len(cube)})  state=2048 bits")
    log(f"active cube bits: {cube[:32]}{' ...' if len(cube)>32 else ''}")
    log("generating S-box division trails ...")
    trails = all_trails()
    log(f"S-box trails ready ({sum(len(v) for v in trails.values())} trails)")
    # Default = exact SELECTOR S-box model (validated, sound). Facets are NOT
    # used: exact facet enumeration of this S-box is explosive (it froze the
    # machine). --use-facets only loads a pre-built facets.pkl; it never
    # auto-generates one.
    facets = None
    if args.use_facets:
        import pickle, os
        pkl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "facets.pkl")
        if not os.path.exists(pkl):
            log("ERROR: --use-facets given but facets.pkl not present (and we never "
                "auto-generate it). Run with the default selector model instead.")
            sys.exit(2)
        with open(pkl, "rb") as f:
            info = pickle.load(f)
        facets = info["facets"]
        log(f"using FACET S-box model: {info['n_facets']} inequalities (loose, sound)")
    else:
        log("using SELECTOR S-box model (exact, validated, ~2015 vars/sbox)")
    log(f"building {args.rounds}-round division MILP (this is the slow part) ...")
    c, out = build(args.rounds, cube, trails, facets=facets)
    nbits = model.NB if args.max_bits <= 0 else min(args.max_bits, model.NB)
    test_bits = range(nbits)
    if args.solver_time_limit > 0:
        c.m.setParam("limits/time", args.solver_time_limit)
        log(f"per-bit SCIP time limit = {args.solver_time_limit}s")
    log(f"model built: {c.m.getNVars()} vars, {c.m.getNConss()} constraints")
    log(f"testing {nbits} output bits for provable balance ...")
    balanced, unknown = [], []
    t0 = time.time()
    for n, j in enumerate(test_bits):
        cons = [c.m.addCons(out[i] == (1 if i == j else 0)) for i in range(model.NB)]
        ts = time.time()
        c.m.optimize()
        dt = time.time() - ts
        st = c.m.getStatus()
        if st == "infeasible":
            balanced.append(j); tag = "BALANCED"
        elif st == "optimal":
            tag = "reachable"     # not provably balanced
        else:
            unknown.append(j); tag = f"unknown({st})"  # timelimit -> conservatively NOT balanced
        c.m.freeTransform()       # return to problem stage so cons can be removed
        for cc in cons:
            c.m.delCons(cc)
        elapsed = time.time() - t0
        rate = (n + 1) / elapsed
        eta = (nbits - n - 1) / rate if rate > 0 else 0
        log(f"  bit {j:4d} [{n+1:4d}/{nbits}] {tag:14s} solve={dt:6.2f}s "
            f"| balanced so far={len(balanced):4d} | elapsed={elapsed:7.1f}s eta={eta:7.1f}s")
    print(f"# tested {nbits} output bits")
    print(f"BALANCED (provable zero-sum): {len(balanced)} bits -> {balanced[:64]}{' ...' if len(balanced)>64 else ''}")
    if unknown:
        print(f"UNKNOWN (hit time limit): {len(unknown)} bits")
    if args.verify_empirical:
        verify_empirical(args.rounds, cube, balanced, nbits)
    return balanced

def verify_empirical(rounds, cube, balanced, nbits):
    """Necessary-condition check: XOR-sum the REAL permutation over the cube and
    confirm every tool-flagged balanced bit is empirically zero."""
    d = len(cube)
    if d > 22:
        print(f"# (skip empirical: cube dim {d} too large for 2^d sweep)")
        return
    perm = L._lib.krakken_permute_scalar_rounds
    perm.restype = None
    acc = [0]*model.NB
    for mask in range(1 << d):
        st = L.StateT()
        for bidx, g in enumerate(cube):
            if (mask >> bidx) & 1:
                st[g >> 6] |= (1 << (g & 63))
        perm(st, rounds)
        for w in range(32):
            v = st[w]
            base = w*64
            while v:
                lsb = v & -v; acc[base + lsb.bit_length()-1] ^= 1; v ^= lsb
    emp_balanced = {j for j in range(nbits) if acc[j] == 0}
    bad = [j for j in balanced if j not in emp_balanced]
    print(f"# empirical XOR-sum over 2^{d} cube values:")
    print(f"#   empirically-zero bits among tested: {len(emp_balanced)}/{nbits}")
    print(f"#   tool-balanced bits that are NOT empirically zero: {len(bad)} (MUST be 0)  {bad[:16]}")
    print("#   EMPIRICAL CROSS-CHECK PASSED" if not bad else "#   EMPIRICAL CROSS-CHECK FAILED -- UNSOUND")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, required=True)
    ap.add_argument("--cube", default="word:0", help="word:W | byte:W:B | bits:i,j,.. | dim:D[:seed]")
    ap.add_argument("--max-bits", type=int, default=0, help="test only first N output bits (0=all 2048)")
    ap.add_argument("--solver-time-limit", type=float, default=0.0, help="per-bit SCIP time cap (0=none)")
    ap.add_argument("--verify-empirical", action="store_true", help="cross-check vs real permutation (small cubes)")
    ap.add_argument("--use-facets", action="store_true", help="use pre-built facets.pkl (not generated by default; selector model is the default)")
    run(ap.parse_args())
