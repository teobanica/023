# Unitary moments and order-preserving maps into noncrossing partitions.
# python3 orders.py
# For 0 in D, the moment of (Z^k Z*^k)^p counts the pairings of the points
# (o^k x^k)^p in which every string joins an up-step and a down-step of the same
# edge, and 2 strings cross only at height differences outside D. We count them
# directly and compare with the number of order-preserving maps from the heights
# {0,...,k-1}, ordered by i <= j when j-i is in D, to the lattice NC(p).
from itertools import permutations, product

def crosses(a, b):
    (p, q), (r, s) = sorted(a), sorted(b)
    return p < r < q < s or r < p < s < q

def direct_count(D, k, p):
    ups = {h: [] for h in range(k)}
    downs = {h: [] for h in range(k)}
    pos = 0
    for _ in range(p):
        for h in range(k):
            ups[h].append(pos); pos += 1
        for h in reversed(range(k)):
            downs[h].append(pos); pos += 1
    choices = []
    for h in range(k):
        level = []
        for perm in permutations(range(p)):
            strings = [(ups[h][i], downs[h][perm[i]]) for i in range(p)]
            if 0 in D and any(crosses(x, y) for x in strings for y in strings if x < y):
                continue
            level.append(strings)
        choices.append(level)
    total = 0
    for pick in product(*choices):
        ok = True
        for i in range(k):
            for j in range(i+1, k):
                if (j-i) in D and any(crosses(x, y) for x in pick[i] for y in pick[j]):
                    ok = False; break
            if not ok: break
        total += ok
    return total

def noncrossing_partitions(p):
    def partitions(items):
        if not items:
            yield []; return
        first, rest = items[0], items[1:]
        for part in partitions(rest):
            yield [[first]] + part
            for i in range(len(part)):
                yield part[:i] + [[first] + part[i]] + part[i+1:]
    out = []
    for part in partitions(list(range(p))):
        blocks = [sorted(b) for b in part]
        bad = False
        for b in blocks:
            for c in blocks:
                if b is not c and any(x < y < z < w for x in b for z in b for y in c for w in c):
                    bad = True
        if not bad:
            out.append(frozenset(frozenset(b) for b in blocks))
    return out

def finer(s, t):
    return all(any(b <= c for c in t) for b in s)

def order_count(D, k, p):
    nc = noncrossing_partitions(p)
    total = 0
    for f in product(range(len(nc)), repeat=k):
        if all(finer(nc[f[i]], nc[f[j]]) for i in range(k) for j in range(i+1, k) if (j-i) in D):
            total += 1
    return total

def semigroup(gens, top=40):
    S = {0}
    changed = True
    while changed:
        changed = False
        for s in list(S):
            for g in gens:
                if s+g <= top and s+g not in S:
                    S.add(s+g); changed = True
    return S

examples = {"N0": semigroup([1]), "{0}": {0}, "D* = <2,3>": semigroup([2, 3]),
            "2N0": semigroup([2]), "<2,5>": semigroup([2, 5]), "<3,4>": semigroup([3, 4]),
            "<3,5,7>": semigroup([3, 5, 7])}
print("D, k, p, pairings, order-preserving maps")
for name, D in examples.items():
    for k, p in [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (1, 3), (2, 3), (3, 3), (4, 3), (2, 4), (3, 4)]:
        a, b = direct_count(D, k, p), order_count(D, k, p)
        print(name, k, p, a, b, "" if a == b else "MISMATCH")
