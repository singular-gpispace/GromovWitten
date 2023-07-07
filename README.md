# TropicalFeynman

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/singular-gpispace/tropicalfeynman/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/singular-gpispace/tropicalfeynman/dev/)
[![Build Status](https://github.com/singular-gpispace/tropicalfeynman/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/singular-gpispace/tropicalfeynman/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/singular-gpispace/tropicalfeynman.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/singular-gpispace/tropicalfeynman.jl)

# TropicalFeynman

The package tropicalfeynman computes generating series for tropical Hurwitz numbers of elliptic curves via mirror symmetry and Feynman integrals, and thus, by a correspondence theorem, Hurwitz numbers in the sense algebraic geometry. Generalizations of the method also allow for the computation of Gromov-Witten invariants for ellptic curves, and are also implemented in the package. TropicalFeynman is based on the computeralgebra system OSCAR and is provided as a package for the Julia programming language.  

# Installation

We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:

```bash
git pull 
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
using TropicalFeynman  
```

# Example

To provide an example on how to use our package, we efine a graph G from a list of edges:

```bash
ve = [ (1, 2), (2, 3), (3, 1)]  
G = graph(ve)
```

We then define a polynomial ring with all variables required by our implementation:

```bash
R,x,q,z=polynomial_ring(G,"x","q","z")
```

Here, the indexed variables x correspond to the vertices of the graph, the indexed variables y to the edges of the graph, and the indexed variables z again to the vertices of the graph (the latter to be used in the context of Gromov-Witten invariants with non-trivial Psi-classes).

To compute a Feynman iuntegral, we define a partition  $a=[0,0,3]$  of degree d=3, a fixed order of vertex $o=[1,2,3]$ and the genus function $g=[1,0,0]$. The leak in G is $L=[0,0,0]$ , $aa=1$ is the order of the Sfunction. We have then

```bash
 feynman_integral_branchtype_order(R,x,q,z,G,a,o,aa=1,l=[0,0,0],g=[1,0,0])
```

also we can compute Feynman Integral of degree 4

```bash
 feynman_integral_degree(R,x,q,z,G,3,aa=1,l=[0,0,0],g=[1,0,0])
```

Finally we substitute all $q$  variables by $q_{1}$

```bash
substitute(feynman_integral_degree(R,x,q,z,G,3,aa=1,l=[0,0,0],g=[1,0,0]))
```
