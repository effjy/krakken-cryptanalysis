"""Gate the AND-rule + carry-chain DIVISION semantics for the modular adder.
Two-subset CBDP is an OVER-approximation: it must be SOUND (exact-reachable bit
=> MILP-reachable bit), never claiming a balanced bit that isn't. We verify
soundness per output bit and also report looseness (MILP reachable but exact not).
"""
from cbdp import CBDP
from model import DivEngine

def addN(eng, x, y, n):
    out=[]; carry=eng.const0()
    for i in range(n):
        xs=eng.copy(x[i],3); ys=eng.copy(y[i],3); cs=eng.copy(carry,3)
        s=eng.xor([xs[0],ys[0],cs[0]])
        m1=eng.and_([xs[1],ys[1]]); m2=eng.and_([xs[2],cs[1]]); m3=eng.and_([ys[2],cs[2]])
        carry=eng.xor([m1,m2,m3]); out.append(s)
    return out

def covers(w,u): return (u & ~w)==0

def exact_bit_reach(n):
    """reach[u] = set of output bit indices j with e_j exactly reachable from u."""
    IN=2*n; SZ=1<<IN
    z=[[0]*SZ for _ in range(n)]
    for a in range(SZ):
        x=a&((1<<n)-1); y=a>>n; s=(x+y)&((1<<n)-1)
        for i in range(n): z[i][a]=(s>>i)&1
    def anf(tt):
        f=tt[:]
        for i in range(IN):
            st=1<<i
            for a in range(SZ):
                if a&st: f[a]^=f[a^st]
        return {a for a in range(SZ) if f[a]}
    bitmons=[anf(z[j]) for j in range(n)]
    reach={}
    for u in range(SZ):
        reach[u]={j for j in range(n) if any(covers(w,u) for w in bitmons[j])}
    return reach

def milp_bit_reachable(n,u,j):
    c=CBDP(); eng=DivEngine(c,None)
    x=c.bits(n); y=c.bits(n); zz=addN(eng,x,y,n)
    for i in range(n):
        c.fix(x[i],(u>>i)&1); c.fix(y[i],(u>>(n+i))&1)
    for i in range(n):
        c.fix(zz[i], 1 if i==j else 0)
    return c.feasible()

if __name__=="__main__":
    n=4
    ex=exact_bit_reach(n)
    unsound=0; loose=0; tot=0
    for u in range(1<<(2*n)):
        for j in range(n):
            tot+=1
            milp=milp_bit_reachable(n,u,j)
            exact=(j in ex[u])
            if exact and not milp:
                unsound+=1
                if unsound<=8: print(f"  UNSOUND u={u:#05x} j={j}: exact reachable but MILP says balanced")
            if milp and not exact:
                loose+=1
    print(f"checked {tot} (u,bit) pairs at n={n}")
    print(f"unsound (exact-reachable but MILP-balanced): {unsound}  <-- MUST be 0")
    print(f"loose   (MILP-reachable but exact-balanced): {loose}  (over-approx; sound but not tight)")
    print("ARX-RULE GATE PASSED (sound)" if unsound==0 else "ARX-RULE GATE FAILED (UNSOUND)")
