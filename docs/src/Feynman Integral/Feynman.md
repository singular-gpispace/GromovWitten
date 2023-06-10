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

A Feynman graph Γ of genus g is a trivalent connected graph of genus g.
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
julia> G=graphe(ve)
graphe([(1, 1), (1, 2), (2, 3), (3, 1)])
```

We then define the $R,x,q=polynomialring(G)$ from the graph G.  The polynomial ring has $ 5g-5$ variables, consisting of two sets of variables: $x_{1},x_{2},...,x_{2g-2}$ and $q_{1},q_{2},...,q_{3g-3}$.

```jldoctest graph
julia>   R,x,q,z=polynomialring(G,"x","q","z")
(Multivariate polynomial ring in 10 variables over QQ, QQMPolyRingElem[x_{1}, x_{2}, x_{3}], QQMPolyRingElem[q_{1}, q_{2}, q_{3}, q_{4}], QQMPolyRingElem[z_{1}, z_{2}, z_{3}])

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
julia>  specificFeynmanIntegralo(R,x,q,z,G,a,o,aa=0,l=[0,0,0],g=[0,0,0])
56*q[2]^4*q[3]^2*q[6]^2
```

## Feynman Integral

The Feynman Integral of the graph G, calculated to a specified degree d.

```jldoctest graph
julia> feynmanIntegral(R,x,q,G,4) # here d=4
352*q[1]^8 + 992*q[1]^6*q[2]^2 + 120*q[1]^6*q[5]^2 + 120*q[1]^6*q[6]^2 + 1576*q[1]^4*q[2]^4 + 256*q[1]^4*q[2]^2*q[5]^2 + 256*q[1]^4*q[2]^2*q[6]^2 + 96*q[1]^4*q[3]^4 + 56*q[1]^4*q[3]^2*q[5]^2 + 56*q[1]^4*q[3]^2*q[6]^2 + 96*q[1]^4*q[4]^4 + 56*q[1]^4*q[4]^2*q[5]^2 + 56*q[1]^4*q[4]^2*q[6]^2 + 96*q[1]^4*q[5]^4 + 96*q[1]^4*q[5]^2*q[6]^2 + 96*q[1]^4*q[6]^4 + 992*q[1]^2*q[2]^6 + 256*q[1]^2*q[2]^4*q[5]^2 + 256*q[1]^2*q[2]^4*q[6]^2 + 32*q[1]^2*q[2]^2*q[3]^4 + 32*q[1]^2*q[2]^2*q[4]^4 + 96*q[1]^2*q[2]^2*q[5]^4 + 48*q[1]^2*q[2]^2*q[5]^2*q[6]^2 + 96*q[1]^2*q[2]^2*q[6]^4 + 432*q[1]^2*q[3]^6 + 136*q[1]^2*q[3]^4*q[5]^2 + 136*q[1]^2*q[3]^4*q[6]^2 + 48*q[1]^2*q[3]^2*q[4]^2*q[5]^2 + 48*q[1]^2*q[3]^2*q[4]^2*q[6]^2 + 56*q[1]^2*q[3]^2*q[5]^4 + 56*q[1]^2*q[3]^2*q[6]^4 + 432*q[1]^2*q[4]^6 + 136*q[1]^2*q[4]^4*q[5]^2 + 136*q[1]^2*q[4]^4*q[6]^2 + 56*q[1]^2*q[4]^2*q[5]^4 + 56*q[1]^2*q[4]^2*q[6]^4 + 120*q[1]^2*q[5]^6 + 256*q[1]^2*q[5]^4*q[6]^2 + 256*q[1]^2*q[5]^2*q[6]^4 + 120*q[1]^2*q[6]^6 + 352*q[2]^8 + 120*q[2]^6*q[5]^2 + 120*q[2]^6*q[6]^2 + 96*q[2]^4*q[3]^4 + 56*q[2]^4*q[3]^2*q[5]^2 + 56*q[2]^4*q[3]^2*q[6]^2 + 96*q[2]^4*q[4]^4 + 56*q[2]^4*q[4]^2*q[5]^2 + 56*q[2]^4*q[4]^2*q[6]^2 + 96*q[2]^4*q[5]^4 + 96*q[2]^4*q[5]^2*q[6]^2 + 96*q[2]^4*q[6]^4 + 432*q[2]^2*q[3]^6 + 136*q[2]^2*q[3]^4*q[5]^2 + 136*q[2]^2*q[3]^4*q[6]^2 + 48*q[2]^2*q[3]^2*q[4]^2*q[5]^2 + 48*q[2]^2*q[3]^2*q[4]^2*q[6]^2 + 56*q[2]^2*q[3]^2*q[5]^4 + 56*q[2]^2*q[3]^2*q[6]^4 + 432*q[2]^2*q[4]^6 + 136*q[2]^2*q[4]^4*q[5]^2 + 136*q[2]^2*q[4]^4*q[6]^2 + 56*q[2]^2*q[4]^2*q[5]^4 + 56*q[2]^2*q[4]^2*q[6]^4 + 120*q[2]^2*q[5]^6 + 256*q[2]^2*q[5]^4*q[6]^2 + 256*q[2]^2*q[5]^2*q[6]^4 + 120*q[2]^2*q[6]^6 + 3208*q[3]^8 + 432*q[3]^6*q[5]^2 + 432*q[3]^6*q[6]^2 + 48*q[3]^4*q[4]^4 + 96*q[3]^4*q[5]^4 + 32*q[3]^4*q[5]^2*q[6]^2 + 96*q[3]^4*q[6]^4 + 3208*q[4]^8 + 432*q[4]^6*q[5]^2 + 432*q[4]^6*q[6]^2 + 96*q[4]^4*q[5]^4 + 32*q[4]^4*q[5]^2*q[6]^2 + 96*q[4]^4*q[6]^4 + 352*q[5]^8 + 992*q[5]^6*q[6]^2 + 1576*q[5]^4*q[6]^4 + 992*q[5]^2*q[6]^6 + 352*q[6]^8
```

We then compute the sum of the feynman integral from the degree 1 to d

```jldoctest graph
julia> feynmanIntegralSum(R,x,q,G,4)
352*q[1]^8 + 992*q[1]^6*q[2]^2 + 120*q[1]^6*q[5]^2 + 120*q[1]^6*q[6]^2 + 24*q[1]^6 + 1576*q[1]^4*q[2]^4 + 256*q[1]^4*q[2]^2*q[5]^2 + 256*q[1]^4*q[2]^2*q[6]^2 + 152*q[1]^4*q[2]^2 + 96*q[1]^4*q[3]^4 + 56*q[1]^4*q[3]^2*q[5]^2 + 56*q[1]^4*q[3]^2*q[6]^2 + 96*q[1]^4*q[4]^4 + 56*q[1]^4*q[4]^2*q[5]^2 + 56*q[1]^4*q[4]^2*q[6]^2 + 96*q[1]^4*q[5]^4 + 96*q[1]^4*q[5]^2*q[6]^2 + 8*q[1]^4*q[5]^2 + 96*q[1]^4*q[6]^4 + 8*q[1]^4*q[6]^2 + 992*q[1]^2*q[2]^6 + 256*q[1]^2*q[2]^4*q[5]^2 + 256*q[1]^2*q[2]^4*q[6]^2 + 152*q[1]^2*q[2]^4 + 32*q[1]^2*q[2]^2*q[3]^4 + 32*q[1]^2*q[2]^2*q[4]^4 + 96*q[1]^2*q[2]^2*q[5]^4 + 48*q[1]^2*q[2]^2*q[5]^2*q[6]^2 + 32*q[1]^2*q[2]^2*q[5]^2 + 96*q[1]^2*q[2]^2*q[6]^4 + 32*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[2]^2 + 432*q[1]^2*q[3]^6 + 136*q[1]^2*q[3]^4*q[5]^2 + 136*q[1]^2*q[3]^4*q[6]^2 + 32*q[1]^2*q[3]^4 + 48*q[1]^2*q[3]^2*q[4]^2*q[5]^2 + 48*q[1]^2*q[3]^2*q[4]^2*q[6]^2 + 56*q[1]^2*q[3]^2*q[5]^4 + 8*q[1]^2*q[3]^2*q[5]^2 + 56*q[1]^2*q[3]^2*q[6]^4 + 8*q[1]^2*q[3]^2*q[6]^2 + 432*q[1]^2*q[4]^6 + 136*q[1]^2*q[4]^4*q[5]^2 + 136*q[1]^2*q[4]^4*q[6]^2 + 32*q[1]^2*q[4]^4 + 56*q[1]^2*q[4]^2*q[5]^4 + 8*q[1]^2*q[4]^2*q[5]^2 + 56*q[1]^2*q[4]^2*q[6]^4 + 8*q[1]^2*q[4]^2*q[6]^2 + 120*q[1]^2*q[5]^6 + 256*q[1]^2*q[5]^4*q[6]^2 + 8*q[1]^2*q[5]^4 + 256*q[1]^2*q[5]^2*q[6]^4 + 32*q[1]^2*q[5]^2*q[6]^2 + 120*q[1]^2*q[6]^6 + 8*q[1]^2*q[6]^4 + 352*q[2]^8 + 120*q[2]^6*q[5]^2 + 120*q[2]^6*q[6]^2 + 24*q[2]^6 + 96*q[2]^4*q[3]^4 + 56*q[2]^4*q[3]^2*q[5]^2 + 56*q[2]^4*q[3]^2*q[6]^2 + 96*q[2]^4*q[4]^4 + 56*q[2]^4*q[4]^2*q[5]^2 + 56*q[2]^4*q[4]^2*q[6]^2 + 96*q[2]^4*q[5]^4 + 96*q[2]^4*q[5]^2*q[6]^2 + 8*q[2]^4*q[5]^2 + 96*q[2]^4*q[6]^4 + 8*q[2]^4*q[6]^2 + 432*q[2]^2*q[3]^6 + 136*q[2]^2*q[3]^4*q[5]^2 + 136*q[2]^2*q[3]^4*q[6]^2 + 32*q[2]^2*q[3]^4 + 48*q[2]^2*q[3]^2*q[4]^2*q[5]^2 + 48*q[2]^2*q[3]^2*q[4]^2*q[6]^2 + 56*q[2]^2*q[3]^2*q[5]^4 + 8*q[2]^2*q[3]^2*q[5]^2 + 56*q[2]^2*q[3]^2*q[6]^4 + 8*q[2]^2*q[3]^2*q[6]^2 + 432*q[2]^2*q[4]^6 + 136*q[2]^2*q[4]^4*q[5]^2 + 136*q[2]^2*q[4]^4*q[6]^2 + 32*q[2]^2*q[4]^4 + 56*q[2]^2*q[4]^2*q[5]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 56*q[2]^2*q[4]^2*q[6]^4 + 8*q[2]^2*q[4]^2*q[6]^2 + 120*q[2]^2*q[5]^6 + 256*q[2]^2*q[5]^4*q[6]^2 + 8*q[2]^2*q[5]^4 + 256*q[2]^2*q[5]^2*q[6]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 120*q[2]^2*q[6]^6 + 8*q[2]^2*q[6]^4 + 3208*q[3]^8 + 432*q[3]^6*q[5]^2 + 432*q[3]^6*q[6]^2 + 288*q[3]^6 + 48*q[3]^4*q[4]^4 + 96*q[3]^4*q[5]^4 + 32*q[3]^4*q[5]^2*q[6]^2 + 32*q[3]^4*q[5]^2 + 96*q[3]^4*q[6]^4 + 32*q[3]^4*q[6]^2 + 8*q[3]^4 + 3208*q[4]^8 + 432*q[4]^6*q[5]^2 + 432*q[4]^6*q[6]^2 + 288*q[4]^6 + 96*q[4]^4*q[5]^4 + 32*q[4]^4*q[5]^2*q[6]^2 + 32*q[4]^4*q[5]^2 + 96*q[4]^4*q[6]^4 + 32*q[4]^4*q[6]^2 + 8*q[4]^4 + 352*q[5]^8 + 992*q[5]^6*q[6]^2 + 24*q[5]^6 + 1576*q[5]^4*q[6]^4 + 152*q[5]^4*q[6]^2 + 992*q[5]^2*q[6]^6 + 152*q[5]^2*q[6]^4 + 8*q[5]^2*q[6]^2 + 352*q[6]^8 + 24*q[6]^6

```jldoctest graph
julia> subt(R,x,q,feynmanIntegralSum(R,x,q,G,4))

25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4

```
