# Test of the edge bound for difference graphs of numerical semigroups.
# python3 isoperimetry.py 9 18
# For each numerical semigroup with all gaps <= g, and each subset Y of {0,...,L}
# containing 0, we count the pairs of Y whose difference is a gap, and compare with
# the count for the interval {0,...,n-1}, n = |Y|. The conjecture in the section
# "Separation of reflection partners" says that the interval has the largest count.
import sys

def gap_sets(g):
    out = []
    for mask in range(1 << g):
        B = [i+1 for i in range(g) if mask >> i & 1]
        gaps = set(B)
        if all(s in gaps or (b-s) in gaps for b in B for s in range(1, b)):
            out.append(B)
    return out

g = int(sys.argv[1]) if len(sys.argv) > 1 else 6
L = int(sys.argv[2]) if len(sys.argv) > 2 else 12
tested = failures = extremal = 0
for B in gap_sets(g):
    if not B: continue
    gapmask = sum(1 << b for b in B)
    interval = [sum(n-b for b in B if b < n) for n in range(L+3)]     # pairs of {0..n-1}
    for half in range(1 << L):
        Y = 2*half + 1                       # bit i set when i belongs to Y; 0 always does
        n = bin(Y).count("1")
        edges = 0
        for b in B:
            edges += bin(Y & (Y >> b)).count("1")
        tested += 1
        if edges > interval[n]:
            failures += 1; print("edge bound fails", B, [i for i in range(L+1) if Y >> i & 1])
        elif edges == interval[n]:
            extremal += 1
print("gaps up to", g, "subsets of 0..%d:" % L, tested, "sets tested,", failures,
      "failures,", extremal, "sets with as many pairs as the interval")
