// Compute the first difference of the limiting Haar moments for the path partner.
// c++ -O2 -std=c++20 path-moments.cpp -o path-moments && ./path-moments
// All word tests are exact. The large-alphabet search is for a shortest new
// relation: deleting a letter must give the identity freely, and each letter
// occurs at least four times. See the proof in "Reflection character laws".
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using Word = std::string;
using Mask = unsigned;
using Signature = std::array<std::uint64_t, 3>;
constexpr int cutoff = 30;

int first_bit(Mask mask) { return std::countr_zero(mask); }

Word rename(const Word& w) {
    std::array<char, 8> names{};
    Word result;
    char next = 'a';
    for (char c : w) {
        char& name = names[c-'a'];
        if (!name) name = next++;
        result += name;
    }
    return result;
}

bool path_identity(const Word& w, const std::vector<int>& site) {
    std::array<int, cutoff> reduced{};
    int size = 0;
    for (char c : w) {
        int v = site[c-'a'];
        int j = size-1;
        // Commute past distinct, nonneighboring sites. A neighbor blocks the
        // cancellation. Equal generators cancel since each is an involution.
        while (j >= 0 && std::abs(reduced[j]-v) > 1) --j;
        if (j >= 0 && reduced[j] == v) {
            std::move(reduced.begin()+j+1, reduced.begin()+size, reduced.begin()+j);
            --size;
        } else {
            reduced[size++] = v;
        }
    }
    return size == 0;
}

bool three_paths(const Word& w) {
    std::vector<int> site = {0, 1, 2};
    do {
        if (!path_identity(w, site)) return false;
    } while (std::next_permutation(site.begin(), site.end()));
    return true;
}

// Base-five digits encode reduced syllables in (Z2 x Z2) * Z2.
// Digits 1,2,3 are the nonzero two-bit masks in the first factor, 4 the middle site.
std::uint64_t append_syllable(std::uint64_t code, int letter, int middle) {
    const int last = code % 5;
    if (letter == middle) return last == 4 ? code/5 : 5*code+4;
    int mask = (letter < middle ? letter : letter-1) == 0 ? 1 : 2;
    if (last == 4 || code == 0) return 5*code+mask;
    mask ^= last;
    return mask == 0 ? code/5 : 5*(code/5)+mask;
}

struct ThreeLetterCheck {
    std::map<Signature, std::vector<Word>> bins;
    std::set<Word> patterns;
    int ordered_pairs = 0;

    void halves(int length, Word& half, const Signature& sig) {
        if (static_cast<int>(half.size()) == length) {
            auto& previous = bins[sig];
            for (const auto& other : previous) {
                Word w = half + Word(other.rbegin(), other.rend());
                patterns.insert(rename(w));
                std::reverse(w.begin(), w.end());
                patterns.insert(rename(w));
                ordered_pairs += 2;
            }
            previous.push_back(half);
            return;
        }
        for (int x = 0; x < 3; ++x) {
            char c = 'a'+x;
            if (!half.empty() && half.back() == c) continue;
            Signature next;
            for (int j = 0; j < 3; ++j) next[j] = append_syllable(sig[j],x,j);
            half.push_back(c);
            halves(length, half, next);
            half.pop_back();
        }
    }
};

class LargerAlphabetCheck {
public:
    LargerAlphabetCheck(int length, int support, int multiplicity, int triple_cutoff)
        : n(length), m(support), l(multiplicity), pairs((n-l)/2),
          first_triple(triple_cutoff), gaps(l), skeleton(n), incompatible(pairs),
          allowed(pairs), colors(pairs), independent(1u << pairs) {}

    void run() { choose_gaps(0,n-l); }
    std::int64_t tested = 0, survivors = 0;

private:
    const int n, m, l, pairs, first_triple;
    std::vector<int> gaps, skeleton;
    std::vector<Mask> incompatible;
    std::vector<std::vector<Mask>> allowed;
    std::vector<char> colors;
    std::vector<char> independent;

    bool gap_representative() const {
        // Every cyclic rotation or reversal can be tested using the same placements.
        for (int start = 0; start < l; ++start) for (int step : {1,-1}) {
            for (int j = 0; j < l; ++j) {
                int g = gaps[(start+step*j+2*l)%l];
                if (g != gaps[j]) {
                    if (g > gaps[j]) return false;
                    break;
                }
            }
        }
        return true;
    }

    void choose_gaps(int j, int remaining) {
        if (j == l) {
            if (!remaining && gap_representative()) choose_walk(0,0,0,{});
            return;
        }
        int limit = remaining-(l-j-1);
        if (j) limit = std::min(limit, gaps[0]);
        for (int g = 1; g <= limit; ++g) {
            gaps[j] = g;
            choose_gaps(j+1, remaining-g);
        }
    }

    void choose_walk(int j, int pos, int created, const std::vector<int>& stack) {
        if (j == l) {
            if (stack.empty() && created == pairs) color_skeleton();
            return;
        }
        skeleton[pos++] = -1; // the distinguished, least frequent letter
        int max_down = std::min(gaps[j],static_cast<int>(stack.size()));
        for (int down = 0; down <= max_down; ++down) {
            int up = gaps[j]-down;
            if (created+up > pairs) continue;
            if (j == l-1 && (up || down != static_cast<int>(stack.size()))) continue;
            auto next = stack;
            for (int i = 0; i < down; ++i) {
                skeleton[pos+i] = next.back();
                next.pop_back();
            }
            for (int i = 0; i < up; ++i) {
                skeleton[pos+down+i] = created+i;
                next.push_back(created+i);
            }
            choose_walk(j+1,pos+gaps[j],created+up,next);
        }
    }

    void distinct(int a, int b) {
        if (a == b) return;
        incompatible[a] |= 1u<<b;
        incompatible[b] |= 1u<<a;
    }

    bool three_color_test(Mask part) const {
        Word reduced;
        for (int slot : skeleton) {
            char c = slot < 0 ? 'a' : ((part >> slot) & 1u ? 'b' : 'c');
            if (!reduced.empty() && reduced.back() == c) reduced.pop_back();
            else reduced.push_back(c);
        }
        if (reduced.empty()) return true;
        return static_cast<int>(reduced.size()) >= first_triple && three_paths(reduced);
    }

    void color_skeleton() {
        std::fill(incompatible.begin(),incompatible.end(),0);
        std::vector<int> stack;
        for (int i = 0; i < n; ++i) {
            int a = skeleton[i];
            if (a < 0) continue;
            if (i && skeleton[i-1] >= 0) distinct(a,skeleton[i-1]);
            if (!stack.empty() && stack.back() == a) stack.pop_back();
            else {
                if (!stack.empty()) distinct(a,stack.back());
                stack.push_back(a);
            }
        }
        for (auto& masks : allowed) masks.clear();
        const Mask all = (1u<<pairs)-1;
        const int largest_part = pairs-(m-2)*(l/2);
        independent[0] = true;
        for (Mask v = 1; v <= all; ++v) {
            int bit = first_bit(v);
            Mask rest = v&(v-1);
            independent[v] = independent[rest] && !(incompatible[bit]&rest);
            int size = std::popcount(v);
            if (independent[v] && size >= l/2 && size <= largest_part && three_color_test(v))
                allowed[bit].push_back(v);
        }
        partition(all,0);
    }

    void partition(Mask remaining, int assigned) {
        if (!remaining) {
            if (assigned != m-1) return;
            ++tested;
            Word w;
            for (int slot : skeleton) w += slot < 0 ? 'a' : colors[slot];
            std::vector<int> site(m);
            std::iota(site.begin(),site.end(),0);
            do {
                if (!path_identity(w,site)) return;
            } while (std::next_permutation(site.begin(),site.end()));
            ++survivors;
            return;
        }
        if (assigned == m-1) return;
        int later = m-2-assigned;
        for (Mask part : allowed[first_bit(remaining)]) {
            if ((remaining&part) != part) continue;
            Mask rest = remaining^part;
            if (std::popcount(rest) < later*(l/2) || (!later && rest)) continue;
            for (Mask bits = part; bits; bits &= bits-1) colors[first_bit(bits)] = 'b'+assigned;
            partition(rest,assigned+1);
        }
    }
};

// Sum m^r (m-1)^(j-r) over Dyck paths, then use the Poisson factorial moments.
// Small integer tables store binomials and Bell numbers through order 15.
std::int64_t tree_moment(int j) {
    std::int64_t binom[31][31]{}, bell[16]{1};
    for (int i = 0; i <= 30; ++i) {
        binom[i][0] = binom[i][i] = 1;
        for (int s = 1; s < i; ++s) binom[i][s] = binom[i-1][s]+binom[i-1][s-1];
    }
    for (int i = 1; i <= j; ++i) for (int s = 0; s < i; ++s)
        bell[i] += binom[i-1][s]*bell[s];
    std::int64_t sum = 0;
    for (int r = 1; r <= j; ++r) {
        std::int64_t weight = 0;
        for (int s = 0; s <= j-r; ++s) {
            auto term = binom[j-r][s]*bell[r+s];
            weight += (j-r-s)%2 ? -term : term;
        }
        sum += r*binom[2*j-r][j]/(2*j-r)*weight;
    }
    return sum;
}

int main() {
    int first_triple = 0, difference = 0;
    for (int half = 1; half <= cutoff/2; ++half) {
        ThreeLetterCheck check;
        Word w;
        check.halves(half,w,{});
        if (check.ordered_pairs) {
            first_triple = 2*half;
            difference = static_cast<int>(check.patterns.size());
            std::cout << "Three letters: length " << first_triple << ", labeled words "
                      << check.ordered_pairs << ", equality patterns " << difference << '\n';
            Word seed = "abacababcbcacbcacababcbabcbcac";
            std::set<Word> rotations;
            for (int i = 0; i < first_triple; ++i) {
                rotations.insert(rename(seed));
                std::rotate(seed.begin(),seed.begin()+1,seed.end());
            }
            if (rotations != check.patterns) throw std::runtime_error("rotation check failed");
            break;
        }
    }
    if (first_triple != cutoff) throw std::runtime_error("unexpected three-letter cutoff");
    std::cout << "Length, candidate colorings on 4,5,6,7 letters, survivors\n";
    for (int n = 16; n <= first_triple; n += 2) {
        std::cout << n;
        std::int64_t survivors = 0;
        for (int m = 4; m <= 7; ++m) {
            std::int64_t tested = 0;
            for (int l = 4; l <= n/m; l += 2) {
                LargerAlphabetCheck check(n,m,l,first_triple);
                check.run();
                tested += check.tested;
                survivors += check.survivors;
            }
            std::cout << " " << tested;
        }
        std::cout << " " << survivors << '\n';
        if (survivors) throw std::runtime_error("check the surviving candidates");
    }
    for (int j = 1; j <= cutoff/2; ++j) {
        auto free = tree_moment(j);
        std::cout << "degree " << 2*j << ": " << free << " "
                  << free + (2*j == first_triple ? difference : 0) << '\n';
    }
}
