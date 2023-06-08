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
Once julia opened type:
```bash

 pkg> instantiate 
```

```bash
using tropicalfeynman  #this will load the package 
```
then we define a graphe using a list.

```bash
 ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
 #this is the Cartepillar graph.  
  G=graphe(ve) #The graphe G.

```
We then define the Polynomial Ring  

```bash
R,x,q,z=polynomialring(G) 
 x represents a vector of vertices,  y represents a vector of edges and  z represents a vector of vertex contributions.
```
To compute the specific Feynman Integral, we define a list ```a=[0,2,1,0,0,1]```;
We have then 
```bash
 specificFeynmanIntegral(R,x,q,G,a)
```
also we can compute Feynman Integral of degree 4
```bash
feynmanIntegral(R,x,q,G ,4)
```
Finally we substitute all q variable by

```bash
subtV(feynmanIntegral(R,x,q,G ,4))
```
