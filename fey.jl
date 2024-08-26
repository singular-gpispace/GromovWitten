using GromovWitten
ve = [(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)]
ve = [(1, 2), (1, 2), (1, 3), (2, 4), (3, 4), (3, 5), (4, 6), (5, 6), (5, 6)]
ve = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
G = FeynmanGraph(ve)
F = FeynmanIntegral(G)
d = 8
@time substitute(feynman_integral_degree(F, d))

