Programs for the paper "Crossing distances and quantum reflection groups".

path-moments.cpp   Section "Reflection character laws". First difference between the
                   limiting moments of the path partner and of the free reflection group.
                   c++ -O2 -std=c++20 path-moments.cpp -o path-moments && ./path-moments

dimension.cpp      Section "Separation of reflection partners". For all pairs of numerical
                   semigroups with gaps at most g, checks that the graph on the heights
                   0,...,k coming from one semigroup does not fit in the graph of the other.
                   c++ -O2 -std=c++17 dimension.cpp -o dimension && ./dimension 12
                   (about 1 minute up to g=11, and several minutes more for g=12)

isoperimetry.py    Section "Separation of reflection partners". Test of the edge bound
                   conjecture. python3 isoperimetry.py 9 18   (about 2 minutes)

orders.py          Section "Unitary character laws". Counts the pairings for the moments of
                   (Z^k Z*^k)^p, and the order-preserving maps into noncrossing partitions.
                   python3 orders.py

cumulants.py       Section "Reflection character laws". Moments and free cumulants of the
                   reflection characters, for several parameters.
                   python3 cumulants.py 8   (order 10 takes about 3 minutes)
