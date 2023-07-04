# tropicalfeynman

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/singular-gpispace/tropicalfeynman/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/singular-gpispace/tropicalfeynman/dev/)
[![Build Status](https://github.com/singular-gpispace/tropicalfeynman/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/singular-gpispace/tropicalfeynman/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/singular-gpispace/tropicalfeynman.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/singular-gpispace/tropicalfeynman.jl)

# tropicalfeynman

To run:

- First, pull the package onto your local disk. Once the package has been successfully pulled, navigate to the desired folder where the package is located, and then type the following command in the terminal..

```bash
julia --project. #this will activate the environment 
```

Once julia opened type this command to install missing packages:

```bash
 pkg> instantiate 
```


or 

```bash
 julia>  import Pkg; Pkg.instantiate()
```


```bash
using tropicalfeynman  #this will load the package 
```

then we define a graph using a  list

```bash
ve=[ (1, 2), (2, 3), (3, 1)]  
G=graph(ve) #The graph G.

```

We then define the Polynomial Ring

```bash
R,x,q,z=polynomial_ring(G,"x","q","z")
 # x represents a vector of vertices, 
 # y represents a vector of edges and 
# z represents a vector of vertex contributions.
```

To compute the specific Feynman Integral, we define a partition  $a=[0,0,3]$  of degree d=3, a fixed order of vertex $o=[1,2,3]$ and the genus function $g=[1,0,0]$. The leak in G is $L=[0,0,0]$ , $aa=1$ is the order of the Sfunction. We have then

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
