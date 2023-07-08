# GromovWitten

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/singular-gpispace/GromovWitten/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/singular-gpispace/GromovWitten/dev/)
[![Build Status](https://github.com/singular-gpispace/GromovWitten/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/singular-gpispace/GromovWitten/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/singular-gpispace/GromovWitten.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/singular-gpispace/GromovWitten.jl)

# GromovWitten

The package GromovWitten computes generating series for tropical Hurwitz numbers of elliptic curves via mirror symmetry and Feynman integrals, and thus, by a correspondence theorem, Hurwitz numbers in the sense algebraic geometry. Generalizations of the method also allow for the computation of Gromov-Witten invariants for ellptic curves, and are also implemented in the package. GromovWitten is based on the computeralgebra system OSCAR and is provided as a package for the Julia programming language.

# Installation

We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:

```bash
git pull https://github.com/singular-gpispace/GromovWitten.git
```

I the same folder execute the following command:

```bash
julia --project
```

This will activate the environment for our package. In Julia install missing packages:

```bash
import Pkg; Pkg.instantiate()
```

and load our package. On the first run this may take some time.

```bash
using GromovWitten  
```

# Example of graph without vertex contribution and loop.

<img width="400" alt="image" src="https://github.com/singular-gpispace/GromovWitten/assets/46294807/0b1f5684-3550-41ea-9722-8403cd96ed35">

```bash
julia> G = graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
```
```bash
julia> R, x, q = polynomial_ring(G, "x", "q")
(Multivariate polynomial ring in 10 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3],x[4]], QQMPolyRingElem[q[1], q[2],q[3], q[4], q[5], q[6]])
```
```bash
julia> a = [0, 2, 1, 0, 0, 1];
```
```bash
julia> o=[1,3,4,2];
```

```bash

julia>  feynman_integral_branchtype_order(x,q,G,a,o) 
128*q[2]^2*q[3]*q[6]
```

```bash
julia> feynman_integral_branchtype(x, q, G, a)  
256*q[2]^2*q[3]*q[6]
```
```bash
julia> f = feynman_integral_degree(x, q, G, 3)
288*q[1]^3 + 32*q[1]^2*q[2] + 32*q[1]^2*q[3] + 32*q[1]^2*q[5] + 32*q[1]^2*q[6] + 8*q[1]*q[2]*q[5] + 8*q[1]*q[2]*q[6] + 8*q[1]*q[3]*q[5] + 8*q[1]*q[3]*q[6] + 24*q[2]^3 + 152*q[2]^2*q[3] + 8*q[2]^2*q[5] + 8*q[2]^2*q[6] + 152*q[2]*q[3]^2 + 32*q[2]*q[3]*q[5] + 32*q[2]*q[3]*q[6] + 32*q[2]*q[4]^2 + 8*q[2]*q[4]*q[5] + 8*q[2]*q[4]*q[6] + 8*q[2]*q[5]^2 + 32*q[2]*q[5]*q[6] + 8*q[2]*q[6]^2 + 24*q[3]^3 + 8*q[3]^2*q[5] + 8*q[3]^2*q[6] + 32*q[3]*q[4]^2 + 8*q[3]*q[4]*q[5] + 8*q[3]*q[4]*q[6] + 8*q[3]*q[5]^2 + 32*q[3]*q[5]*q[6] + 8*q[3]*q[6]^2 + 288*q[4]^3 + 32*q[4]^2*q[5] + 32*q[4]^2*q[6] + 24*q[5]^3 + 152*q[5]^2*q[6] + 152*q[5]*q[6]^2 + 24*q[6]^3
```
```bash
julia>     substitute(q,feynman_integral_degree_sum(x,q,G,8))
10246144*q[1]^8 + 3294720*q[1]^7 + 886656*q[1]^6 + 182272*q[1]^5 + 25344*q[1]^4 + 1792*q[1]^3 + 32*q[1]^2
```
# Example of graph without vertex contribution

<img width="400" alt="image" src="https://github.com/singular-gpispace/GromovWitten/assets/46294807/e5ed2790-64f4-4853-a99c-61b082ddfd73">

To provide an example on how to use our package, we efine a graph G from a list of edges:

```bash
julia> ve = [ (1, 2), (2, 3), (3, 1)]  
julia> G = graph(ve)
```
We then define a polynomial ring with all variables required by our implementation:

```bash
julia> R,x,q,z=polynomial_ring(G,"x","q","z")
```


Here, the indexed variables x correspond to the vertices of the graph, the indexed variables y to the edges of the graph, and the indexed variables z again to the vertices of the graph (the latter to be used in the context of Gromov-Witten invariants with non-trivial Psi-classes).

To compute a Feynman iuntegral, we define a partition  $a=[0,0,3]$  of degree d=3, a fixed order of vertex $o=[1,2,3]$ and the genus function $g=[1,0,0]$. The leak in G is $L=[0,0,0]$ , $aa=1$ is the order of the Sfunction. We have then

```bash
 julia> feynman_integral_branchtype_order(R,x,q,z,G,a,o,aa=1,l=[0,0,0],g=[1,0,0])
```

also we can compute Feynman Integral of degree 4

```bash
 julia> feynman_integral_degree(R,x,q,z,G,3,aa=1,l=[0,0,0],g=[1,0,0])
```

Finally we substitute all $q$  variables by $q_{1}$

```bash
julia> substitute(feynman_integral_degree(R,x,q,z,G,3,aa=1,l=[0,0,0],g=[1,0,0]))
```
# Example of graph with loop.
<img width="350" alt="image" src="https://github.com/singular-gpispace/GromovWitten/assets/46294807/ac17a579-426c-4d16-b652-19cb393d620e">

```bash
julia> G=graph([(1, 1), (1, 2), (2, 3), (3, 1)])
graph([(1, 1), (1, 2), (2, 3), (3, 1)])
```

```bash
julia> R,x,q=polynomial_ring(G,"x","q")
(Multivariate polynomial ring in 7 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3]], QQMPolyRingElem[q[1], q[2], q[3], q[4]])
```

```bash
julia> O=[1,2,3]  
3-element Vector{Int64}:
 1
 2
 3
```

```bash
julia> a=[ 2,  0, 0, 1]
 4-element Vector{Int64}:
 2
 0
 0
 1
```

```bash
julia> feynman_integral_branchtype_order(R,x,q,G,a,O)
3*q[1]^2*q[4]
```

```bash
julia> feynman_integral_branchtype(R,x,q,G,a)  
6*q[1]^2*q[4]
```

```bash
julia> feynman_integral_degree(R,x,q,G,3)
6*q[1]^2*q[2] + 6*q[1]^2*q[3] + 6*q[1]^2*q[4] + 18*q[1]*q[2]^2 + 6*q[1]*q[2]*q[3] + 6*q[1]*q[2]*q[4] + 18*q[1]*q[3]^2 + 6*q[1]*q[3]*q[4] + 18*q[1]*q[4]^2
```

```bash
julia> substitute(R,x,q,feynman_integral_sum(R,x,q,G,8))
20640*q[1]^8 + 9996*q[1]^7 + 4320*q[1]^6 + 1650*q[1]^5 + 456*q[1]^4 + 90*q[1]^3 + 6*q[1]^2
```
