# tropicalfeynman

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/singular-gpispace/tropicalfeynman/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/singular-gpispace/tropicalfeynman/dev/)
[![Build Status](https://github.com/singular-gpispace/tropicalfeynman/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/singular-gpispace/tropicalfeynman/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/singular-gpispace/tropicalfeynman.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/singular-gpispace/tropicalfeynman.jl)

# tropicalfeynman

To run:
- pull the package in your local disk
In the local folder where the package is located, type in the terminal.

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
R,x,q=polynomialring(G) 
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
