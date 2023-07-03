# Feynman Integral

```@meta
CurrentModule = tropicalfeynman
DocTestSetup = quote
  using tropicalfeynman
end
```

```@setup tropicalfeynman
using tropicalfeynman
```

export graphe,Polynomialring
export polynomialring,constterm, proterm, propagator, coefterm, partition, preimg, sgn,flip,flipV

## Graph

A Feynman graph is a (non-metrized) graph Γ without ends with n vertices which are labeled $ x_1, . . . , x_n $ and with labeled edges $q_1, . . . , q_r$.
The graph $G$ is represented as a collection of vertices $V$ and edges $E$. Each edge is a pair $(v,w)$ where both $v$ and $w$ are elements of the set of vertices $V$.
![Bildschirmfoto vom 2023-06-10 17-06-36](https://github.com/singular-gpispace/tropicalfeynman/assets/46294807/c5b4b792-6d2f-418f-b38a-21b3c0187a92)

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
julia> G=graph(ve)
graph([(1, 1), (1, 2), (2, 3), (3, 1)])
```

We then define the $R,x,q=polynomialring(G)$ from the graph G.  The polynomial ring has $ 5g-5$ variables, consisting of two sets of variables: $x_{1},x_{2},...,x_{2g-2}$ and $q_{1},q_{2},...,q_{3g-3}$.

```jldoctest graph
julia>   R,x,q,z=polynomialring(G,"x","q","z")
(Multivariate polynomial ring in 10 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3], q[4]], Nemo.QQMPolyRingElem[z[1], z[2], z[3]])
```

## specific Feynman Integral

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
Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ aa=0$

```jldoctest graph
julia>  specificFeynmanIntegralo(x,q,z,G,a,o,aa=0,l=[0,0,0],g=[0,0,0])
3*q[1]^2*q[4]
```

In the case l=[0,0,0], g=[0,0,0] and  $ aa=0$ we can write simply write.

```jldoctest graph
julia>  specificFeynmanIntegralo(x,q,z,G,a,o)
3*q[1]^2*q[4]
```

We compute the Specific Feynman Integral for the Graph G given a fixed partition of degree $a$ for all vertex ordering $o$.
Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ aa=0$

```jldoctest graph
julia> specificFeynmanIntegral(x,q,z,G,a)
6*q[1]^2*q[4]
```

## Feynman Integral

We compute the  Feynman Integral of the graph G over all  partitions of the degree d=3  for a fixed ordering $o$.

```jldoctest graph
julia> feynmanIntegralo(x,q,z,G,o,3) # here d=3
3*q[1]^2*q[4] + q[1]*q[2]*q[3] + q[1]*q[2]*q[4] + q[1]*q[3]*q[4] + 9*q[1]*q[4]^2

```

We compute the Feynman integral over all the partitions of the degree d of graph G for all vertex ordering.

```jldoctest graph
julia>  feynmanIntegral(x,q,z,G,3) # here d=3
6*q[1]^2*q[2] + 6*q[1]^2*q[3] + 6*q[1]^2*q[4] + 18*q[1]*q[2]^2 + 6*q[1]*q[2]*q[3] + 6*q[1]*q[2]*q[4] + 18*q[1]*q[3]^2 + 6*q[1]*q[3]*q[4] + 18*q[1]*q[4]^2
```

We compute the sum of coefficients.

```
julia> subt( feynmanIntegral(x,q,z,G,3)
90
```
