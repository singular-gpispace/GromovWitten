# Feynman Integral

```@meta
CurrentModule = GromovWitten
DocTestSetup = quote
  using GromovWitten
end
```

```@setup TropicalFeynman
using GromovWitten
```

## Graph

A Feynman graph is a (non-metrized) graph Γ without ends with n vertices which are labeled $x_1, . . . , x_n$ and with labeled edges $q_1, . . . , q_r$.
The graph $G$ is represented as a collection of vertices $V$ and edges $E$. Each edge is a pair $(v,w)$ where both $v$ and $w$ are elements of the set of vertices $V$.

```jldoctest graph
julia> ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
4-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 2)
 (2, 3)
 (3, 1)
```

and

```jldoctest graph
julia> G=FeynmanGraph(ve)
FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)])
```

We define the Feynman integral type, which contains the polynomial Ring and the graph G.

```jldoctest graph
julia> F=FeynmanIntegral(G)
FeynmanIntegral(FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)]), Dict{Symbol, Dict{Vector{Int64}, Nemo.QQMPolyRingElem}}(), (Multivariate polynomial ring in 10 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3], q[4]], Nemo.QQMPolyRingElem[z[1], z[2], z[3]]))
```

## Feynman Integral branche type

For a given branch type $a$, we compute the Specific Feynman Integral of the labeled Graph.
a is a list of partition of degree d=3 of $\Gamma$.

```jldoctest graph
julia> a=[2,0,0,1]
4-element Vector{Int64}:
 2
 0
 0
 1
```

$o$ is a fixed order of vertices.

```jldoctest graph
julia> o=[1,2,3]
3-element Vector{Int64}:
 1
 2
 3
```

We compute the Specific Feynman Integral for the Graph G given a fixed vertex ordering $o$ and the partition of degree $a$.
Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ m=0$

```jldoctest graph
julia>  feynman_integral_branch_type_order(F,a,o)
3*q[1]^4*q[4]^2
```

```We compute the Specific Feynman Integral for the Graph G given a fixed partition of degree $a$ for all vertex ordering $o$.
Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of S-function $ m=0$

```jldoctest graph
julia> feynman_integral_branch_type(F,a)
6*q[1]^4*q[4]^2
```

## Feynman Integral

We compute the  Feynman Integral of the graph G over all  partitions of the degree d=3  for a fixed ordering $o$.

```jldoctest graph
julia> feynman_integral_degree_order(F,o,3) # here d=3
3*q[1]^4*q[4]^2 + q[1]^2*q[2]^2*q[3]^2 + q[1]^2*q[2]^2*q[4]^2 + q[1]^2*q[3]^2*q[4]^2 + 9*q[1]^2*q[4]^4

```

We compute the Feynman integral over all the partitions of the degree d of graph G for all vertex ordering.

```jldoctest graph
julia>  feynman_integral_degree(F,3) # here d=3
6*q[1]^4*q[2]^2 + 6*q[1]^4*q[3]^2 + 6*q[1]^4*q[4]^2 + 18*q[1]^2*q[2]^4 + 6*q[1]^2*q[2]^2*q[3]^2 + 6*q[1]^2*q[2]^2*q[4]^2 + 18*q[1]^2*q[3]^4 + 6*q[1]^2*q[3]^2*q[4]^2 + 18*q[1]^2*q[4]^4
```

We compute the sum of coefficients.

```
julia> sum_of_coeff( feynman_integral_degree(F,3))
90
```
