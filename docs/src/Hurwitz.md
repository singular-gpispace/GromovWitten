# Examples Hurwitz numbers.

![alt text](img/Cartepillar3.png)

To provide an example on how to use our package, we define a graph G from a list of edges:

```julia
julia> G = FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
```

We then define a polynomial ring with all variables required by our implementation:

```julia
julia> F=FeynmanIntegral(G)
```

Here, the indexed variables x correspond to the vertices of the graph, the indexed variables y to the edges of the graph, and the indexed variables z again to the vertices of the graph (the latter to be used in the context of Gromov-Witten invariants with non-trivial Psi-classes).

To compute a Feynman iuntegral, we define a partition  $a = [0, 2, 1, 0, 0, 1]$  of degree d=3, a fixed order of vertex $ o=[1,3,4,2. We have then

```julia
julia> a = [0, 2, 1, 0, 0, 1];
```

```julia
julia> o=[1,3,4,2];
```

The Feynman Integral branch type for a fixed ordering  is

```julia
julia>  feynman_integral_branch_type_order(F,a,o) 
128*q[2]^4*q[3]^2*q[6]^2
```

The Feynman Integral branch type for all ordering  is

```julia
julia> feynman_integral_branch_type(F, a)  
256*q[2]^4*q[3]^2*q[6]^2
```

also we can compute Feynman Integral of degree 3

```julia
julia> f = feynman_integral_degree(F, 3)
288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 24*q[6]^6
```

Finally we substitute all $q$  variables by $q_{1}$ after computing the sum of all Feynman Integral of degree up to 8.

```julia
julia>     substitute(q,feynman_integral_degree_sum(F,8))
10246144*q[1]^16 + 3294720*q[1]^14 + 886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
```
