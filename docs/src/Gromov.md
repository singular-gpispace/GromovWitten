# Examples Gromov-Witten invariants

## Graph with vertex contribution

![alt text](img/graph_with_vertex1.png)

To provide an example on how to use our package, we define a graph G from a list of edges:

```julia
julia> ve = [ (1, 2), (2, 3), (3, 1)]  
julia> G = FeynmanGraph(ve);
```

We then define a polynomial ring with all variables required by our implementation:

```julia
julia> F=FeynmanIntegral(G);
```

To compute a Feynman iuntegral, we define a partition  $a=[0,0,3]$  of degree d=3, a fixed order of vertex $o=[1,2,3]$ and the genus function $g=[1,0,0]$. The leak in G is $L=[0,0,0]$ , $m=1$ is the  default order of the sfunction. We have then

```julia
julia> g=[1,0,0];
julia> a=[0,0,3];
julia> o=[1,2,3]; 
```

```julia
 julia> feynman_integral_branch_type_order(F,a,o,g)
115//6*q[3]^6
```

The Feynman Integral branch type for all ordering with genus function $g$  is

```julia
julia> feynman_integral_branch_type(F,a,g)
115//3*q[3]^6
```

also we can compute Feynman Integral of degree 4

```julia
 julia> feynman_integral_degree(F,4,g)
2041//12*q[1]^8 + 1//4*q[1]^6*q[2]^2 + 1//4*q[1]^6*q[3]^2 + 57//4*q[1]^4*q[2]^4 + 1//2*q[1]^4*q[2]^2*q[3]^2 + 57//4*q[1]^4*q[3]^4 + 1//4*q[1]^2*q[2]^6 + 1//2*q[1]^2*q[2]^4*q[3]^2 + 1//2*q[1]^2*q[2]^2*q[3]^4 + 1//4*q[1]^2*q[3]^6 + 2041//12*q[2]^8 + 1//4*q[2]^6*q[3]^2 + 57//4*q[2]^4*q[3]^4 + 1//4*q[2]^2*q[3]^6 + 2041//12*q[3]^8
```

Finally we substitute all $q$  variables by $q_{1}$

```julia
julia> substitute(feynman_integral_degree(F,4,g))
556*q[1]^8
```

## Graph with loops.

![alt text](img/graph_loop.png)
We have a here a loop at (1,1).

```julia
julia> G=FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)])
FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)])
```

```julia
julia> F=FeynmanIntegral(G);
```

```julia
julia> O=[1,2,3];  # vertices order
```

```julia
julia> a=[ 2,  0, 0, 1] # branch type.
 4-element Vector{Int64}:
 2
 0
 0
 1
```

```julia
julia> feynman_integral_branch_type_order(F, a, O)
3*q[1]^4*q[4]^2
```

```julia
julia> feynman_integral_branch_type(F, a)  
6*q[1]^4*q[4]^2
```

```julia
julia> feynman_integral_degree(F, 3)
6*q[1]^4*q[2]^2 + 6*q[1]^4*q[3]^2 + 6*q[1]^4*q[4]^2 + 18*q[1]^2*q[2]^4 + 6*q[1]^2*q[2]^2*q[3]^2 + 6*q[1]^2*q[2]^2*q[4]^2 + 18*q[1]^2*q[3]^4 + 6*q[1]^2*q[3]^2*q[4]^2 + 18*q[1]^2*q[4]^4
```

```julia
julia> substitute(q,feynman_integral_sum(F, 8))
20640*q[1]^16 + 9996*q[1]^14 + 4320*q[1]^12 + 1650*q[1]^10 + 456*q[1]^8 + 90*q[1]^6 + 6*q[1]^4
```
