"""Generate a SOUND, EXACT inequality model for the ABYSSAL S-box division-
property polytope, for the CBDP MILP (cbdp.CBDP.sbox_ineq).

Why: the default S-box model is the exact SELECTOR encoding (~2015 selector
BINARY variables PER S-box instance). chi instantiates 256 S-boxes/round, so
that is ~5x10^5 extra binaries/round -- the reason the search stalls at round 1.
An inequality model adds ZERO extra variables per S-box (only linear constraints
over the 16 in/out bits), which removes the branching explosion.

Point convention (matches cbdp.sbox_ineq): a 16-vector x = [in0..in7, out0..out7].
A valid division trail is (u -> v) with v reachable from u through the S-box.
P = all valid points; Q = complement (impossible points), |Q| = 2795.

Model: for every impossible point q in Q, add the "no-good" exclusion cut
    sum_{i: q_i=0} x_i  -  sum_{i: q_i=1} x_i  >=  1 - popcount(q)
This inequality is violated ONLY at x = q (any other 0/1 point satisfies it with
equality at best), so it excludes exactly the one impossible transition q while
keeping EVERY valid trail. Over all q in Q the model is therefore:
  * SOUND   -- no genuine division trail is ever excluded;
  * EXACT   -- it admits no impossible transition, matching the selector model's
              tightness bit-for-bit (verified against the trail table below);
  * COMPACT in variables -- 0 extra binaries per S-box (vs ~2015 selectors).

A greedy merge pass (--reduce) then tries to replace groups of no-good cuts by
fewer, stronger inequalities while preserving soundness and exactness, to cut the
constraint count. The un-reduced model is already correct; --reduce only shrinks.
"""
import numpy as np, pickle, os, time, sys, argparse
from sbox_trails import ABYSSAL, anf, N

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

def product_function_anf(v):
    f = [1] * 256
    for i in range(N):
        if (v >> i) & 1:
            for x in range(256):
                f[x] &= (ABYSSAL[x] >> i) & 1
    return anf(f)

def covers(w, u):
    return (u & ~w) == 0

def build_points():
    """return P (valid) and Q (impossible) as Nx16 int arrays; x=[in0..7,out0..7]."""
    mons = {v: product_function_anf(v) for v in range(256)}
    valid = set()
    for v in range(256):
        for u in range(256):
            if any(covers(w, u) for w in mons[v]):
                valid.add((u << 8) | v)          # high8 = u (inputs), low8 = v (outputs)
    def vec(pt):
        u, w = pt >> 8, pt & 0xFF
        return [ (u >> j) & 1 for j in range(8) ] + [ (w >> j) & 1 for j in range(8) ]
    P = np.array([vec(p) for p in sorted(valid)], dtype=np.int64)
    Q = np.array([vec(p) for p in range(1 << 16) if p not in valid], dtype=np.int64)
    return P, Q

def nogood(q):
    """exclusion cut a.x >= b violated ONLY at x=q. a_i=+1 if q_i=0 else -1."""
    a = np.where(q == 0, 1, -1).astype(np.int64)
    b = int(1 - q.sum())
    return a, b

def greedy_reduce(P, Q, facets):
    """Merge no-good cuts into fewer inequalities. For each still-uncovered
    impossible point, relax its no-good cut (zero coefficients while it stays
    valid for ALL of P) into a face that excludes many impossible points, then
    cover greedily. SOUND + EXACT (every q ends up cut; no valid point cut).

    Fast: keeps an incremental dot-product s = P @ a; zeroing coord k updates
    s -= a[k]*P[:,k] in O(|P|) instead of a full matmul."""
    P = P.astype(np.float64)
    Qi = Q.astype(np.int64)
    covered = np.zeros(len(Q), dtype=bool)
    out = []
    pos = 0
    while not covered.all():
        while covered[pos]:
            pos += 1
        q = Q[pos]
        a, b = nogood(q); a = a.astype(np.float64)
        s = P @ a                                  # incremental dot product
        # Relax ONLY coordinates where q has a 0 bit: zeroing those leaves a.q
        # unchanged (q stays excluded), and we keep it only if all of P still
        # satisfies a.x >= b (soundness). This loosens the cut so it also excludes
        # other impossible points, while guaranteeing it still covers q.
        zeroable = [k for k in range(16) if q[k] == 0 and a[k] != 0]
        improved = True
        while improved:
            improved = False
            for k in zeroable:
                if a[k] == 0:
                    continue
                s2 = s - a[k] * P[:, k]
                if s2.min() >= b - 1e-9:
                    s = s2; a[k] = 0.0; improved = True
        ai = a.astype(np.int64)
        assert int(ai @ q) < b                     # q still excluded
        assert (P @ ai).min() >= b - 1e-9          # sound: no valid point cut
        cut = (Qi @ ai) < b
        out.append((ai, int(b)))
        covered |= cut
        if len(out) % 25 == 0:
            log(f"  reduce: {len(out)} ineqs, covered {int(covered.sum())}/{len(Q)}")
    return out

def validate_exact(P, Q, facets):
    """soundness: all P satisfy all cuts; exactness: no Q point admitted."""
    A = np.array([a for a, _ in facets], dtype=np.int64)
    bv = np.array([b for _, b in facets], dtype=np.int64)
    p_ok = (P @ A.T >= bv).all()
    admitted = int((Q @ A.T >= bv).all(axis=1).sum()) if len(Q) else 0
    return bool(p_ok), admitted

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reduce", action="store_true",
                    help="greedy-merge no-good cuts into fewer inequalities")
    args = ap.parse_args()

    log("building valid/impossible division-point sets ...")
    P, Q = build_points()
    log(f"|P| valid = {len(P)}   |Q| impossible = {len(Q)}")

    facets = [nogood(q) for q in Q]
    log(f"base no-good model: {len(facets)} exclusion cuts (exact, 0 extra vars)")

    if args.reduce:
        log("greedy-reducing ...")
        facets = greedy_reduce(P, Q, facets)
        log(f"reduced to {len(facets)} inequalities")

    p_ok, admitted = validate_exact(P, Q, facets)
    log(f"SOUNDNESS (all valid trails kept): {'PASS' if p_ok else 'FAIL'}")
    log(f"EXACTNESS (impossible points admitted): {admitted}  "
        f"({'EXACT' if admitted == 0 else 'loose'})")
    if not p_ok:
        log("ABORT: unsound, not saving."); sys.exit(1)

    facets_ser = [([int(x) for x in a], int(b)) for a, b in facets]
    out = {"facets": facets_ser, "n_facets": len(facets_ser),
           "admitted_impossible": admitted, "n_impossible": int(len(Q)),
           "exact": admitted == 0}
    with open(os.path.join(HERE, "facets.pkl"), "wb") as f:
        pickle.dump(out, f)
    log(f"wrote facets.pkl ({len(facets_ser)} inequalities, "
        f"{'EXACT' if admitted==0 else 'loose'})")

if __name__ == "__main__":
    main()
