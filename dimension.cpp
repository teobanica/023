// First distinguishing dimension for the reflection partners of 2 numerical semigroups.
// c++ -O2 -std=c++17 dimension.cpp -o dimension && ./dimension 12
// For every pair of distinct numerical semigroups S,T with all gaps <= g, let k be the
// smallest integer lying in exactly one of them, say k is a gap of S. The program tests
// whether the graph on {0,...,k}, with i,j joined when |i-j| is a gap of S, can be placed
// injectively in the integers so that joined points have differences which are gaps of T.
// Output: one line for each value g of the largest gap of the pair. See the section
// "Separation of reflection partners".
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

using Gaps = std::vector<int>;

std::vector<Gaps> semigroups(int g) {
    std::vector<Gaps> out;
    for (unsigned mask = 0; mask < (1u << g); ++mask) {
        std::vector<char> gap(2*g+3, 0);
        for (int i = 0; i < g; ++i) if (mask >> i & 1) gap[i+1] = 1;
        bool ok = true;
        for (int b = 1; b <= g && ok; ++b) if (gap[b])
            for (int s = 1; s < b; ++s) if (!gap[s] && !gap[b-s]) { ok = false; break; }
        if (!ok) continue;
        Gaps B;
        for (int b = 1; b <= g; ++b) if (gap[b]) B.push_back(b);
        out.push_back(B);
    }
    return out;
}

struct Embed {
    int k;
    std::vector<char> gapS, gapT;
    std::vector<int> steps, pos;
    std::set<int> used;
    bool rec(int i) {
        if (i > k) return true;
        for (int s : steps) {
            int p = pos[i-1] + s;
            if (used.count(p)) continue;
            bool ok = true;
            for (int j = 0; j < i && ok; ++j)
                if (gapS[i-j]) { int d = std::abs(p - pos[j]); ok = d < (int)gapT.size() && gapT[d]; }
            if (!ok) continue;
            pos[i] = p; used.insert(p);
            if (rec(i+1)) return true;
            used.erase(p);
        }
        return false;
    }
};

int main(int argc, char** argv) {
    int gmax = argc > 1 ? std::atoi(argv[1]) : 10;
    for (int g = 1; g <= gmax; ++g) {
        auto all = semigroups(g);
        std::int64_t pairs = 0, placed = 0; int worst = 0;
        for (const auto& B : all) for (const auto& C : all) {
            if (B == C) continue;
            if ((B.empty() ? 0 : B.back()) != g && (C.empty() ? 0 : C.back()) != g) continue;
            std::vector<char> gb(g+2, 0), gc(g+2, 0);
            for (int b : B) gb[b] = 1;
            for (int c : C) gc[c] = 1;
            int k = 1; while (gb[k] == gc[k]) ++k;
            if (!gb[k]) continue;            // treat the ordered pair where k is a gap of B
            ++pairs; worst = std::max(worst, k+1);
            Embed e; e.k = k; e.gapS = gb; e.gapT = gc; e.pos.assign(k+1, 0); e.used = {0};
            for (int c : C) { e.steps.push_back(c); e.steps.push_back(-c); }
            if (e.rec(1)) { ++placed; std::cout << "placement found, g=" << g << "\n"; }
        }
        std::cout << "g=" << g << ": semigroups " << all.size() << ", pairs " << pairs
                  << ", largest k+1 = " << worst << ", placements found " << placed << std::endl;
    }
}
