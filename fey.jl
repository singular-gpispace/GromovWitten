using GromovWitten
ve = [(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)]
G = FeynmanGraph(ve)
F = FeynmanIntegral(G)
d = 14
@time substitute(feynman_integral_degree(F, d))