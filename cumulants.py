# Reflection moments and free cumulants.
# python3 cumulants.py
# The k-th limiting moment of the main character of the reflection partner counts
# the partitions of {1,...,k} whose word is a relation of the group D'. A word is
# a relation when each letter occurs an even number of times and it reduces to 1
# for every choice of generators s_i of W_D replacing its letters. We count these
# partitions, count the connected ones (no proper noncrossing splitting), and
# check the free moment-cumulant formula.
from itertools import product
import sys

def set_partitions(n):
    def rec(i, blocks):
        if i == n:
            yield [tuple(b) for b in blocks]; return
        for b in blocks:
            b.append(i); yield from rec(i+1, blocks); b.pop()
        blocks.append([i]); yield from rec(i+1, blocks); blocks.pop()
    yield from rec(0, [])

def reduces_to_identity(word, commute):
    stack = []
    for v in word:
        j = len(stack)-1
        while j >= 0 and stack[j] != v and commute(stack[j], v): j -= 1
        if j >= 0 and stack[j] == v: del stack[j]
        else: stack.append(v)
    return not stack

def is_relation(blocks, k, sites, commute):
    if any(len(b) % 2 for b in blocks): return False
    letter = [0]*k
    for a, b in enumerate(blocks):
        for i in b: letter[i] = a
    for f in product(sites, repeat=len(blocks)):
        if 0 not in [s[0] for s in f]: continue          # translation invariance
        if not reduces_to_identity([f[a] for a in letter], commute): return False
    return True

def connected(blocks, k):
    # the noncrossing closure is the 1-block partition
    parent = list(range(len(blocks)))
    def find(x):
        while parent[x] != x: x = parent[x]
        return x
    for a in range(len(blocks)):
        for b in range(a):
            A, B = blocks[a], blocks[b]
            if any(x < y < z < w or y < x < w < z for x in A for z in A for y in B for w in B):
                parent[find(a)] = find(b)
    return len({find(a) for a in range(len(blocks))}) == 1

def noncrossing(k):
    out = []
    for blocks in set_partitions(k):
        if all(not any(x < y < z < w for x in A for z in A for y in B for w in B)
               for A in blocks for B in blocks if A != B):
            out.append([len(b) for b in blocks])
    return out

def study(name, D, twin, window, kmax):
    # sites are pairs (height, copy); twin = 2 when 0 is not in D
    def commute(u, v):
        if u[0] == v[0]: return False                     # 2 copies at one height are free
        return abs(u[0]-v[0]) in D
    moments, cumulants = {}, {}
    for k in range(2, kmax+1, 2):
        sites = [(h, c) for h in range(window*(k//2-1)+1) for c in range(twin)]
        m = c = 0
        for blocks in set_partitions(k):
            if is_relation(blocks, k, sites, commute):
                m += 1; c += connected(blocks, k)
        moments[k], cumulants[k] = m, c
    check = {}
    for k in range(2, kmax+1, 2):
        total = 0
        for sizes in noncrossing(k):
            term = 1
            for s in sizes: term *= cumulants.get(s, 0)
            total += term
        check[k] = total
    print(name, "moments", [moments[k] for k in sorted(moments)],
          "cumulants", [cumulants[k] for k in sorted(cumulants)],
          "formula ok" if check == moments else "MISMATCH %s" % check)

def semigroup(gens, top=60):
    S = {0}
    for _ in range(top):
        S |= {s+g for s in S for g in gens if s+g <= top}
    return S

kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
study("D={0}", {0}, 1, 1, kmax)
study("D=<2,3> (path)", semigroup([2, 3]), 1, 2, kmax)
study("D=<2,5>", semigroup([2, 5]), 1, 4, kmax)
study("D=<3,4,5>", semigroup([3, 4, 5]), 1, 3, kmax)
study("D=2N0", semigroup([2]), 1, 1, kmax)
study("D={2,3,4,...} (doubled path)", semigroup([2, 3]) - {0}, 2, 2, kmax)
study("D=N0", semigroup([1]), 1, 1, kmax)
study("D={1,2,3,...} (single edge)", semigroup([1]) - {0}, 2, 0, kmax)
