```@meta
CurrentModule = GromovWitten
```
# Example of graph with vertex contribution

![alt text](img/graph_with_vertex1.png)

To provide an example on how to use our package, we define a graph G from a list of edges:

```julia
julia> ve = [ (1, 2), (2, 3), (3, 1)]  
julia> G = FeynmanGraph(ve)
```

We then define the Feynman integral type F:

```julia
julia> F=FeynmanIntegral(G)
FeynmanIntegral(FeynmanGraph([(1, 2), (2, 3), (3, 1)]), Dict{Symbol, Dict{Vector{Int64}, Nemo.QQMPolyRingElem}}(), (Multivariate polynomial ring in 9 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3]], Nemo.QQMPolyRingElem[z[1], z[2], z[3]]))
```

Here, the indexed variables x correspond to the vertices of the graph, the indexed variables y to the edges of the graph, and the indexed variables z again to the vertices of the graph (the latter to be used in the context of Gromov-Witten invariants with non-trivial Psi-classes).

To compute a Feynman iuntegral, we define a partition  $a=[0,0,3]$  of degree d=3, a fixed order of vertex $o=[1,2,3]$ and the genus function $g=[1,0,0]$. The leak in G is $L=[0,0,0]$ , $aa=1$ is the order of the sfunction.

```julia
julia> g=[1,0,0];
julia> a=[0,0,3];
julia> o=[1,2,3]; 
```
We have then

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

# Example of graph without vertex contribution and loop.

<img width="400" alt="image" src="https://github.com/singular-gpispace/GromovWitten/assets/46294807/1b45577b-3c92-464f-81f5-57766dcd189e">

```julia
julia> G = FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
```

```julia
julia> F=FeynmanIntegral(G)
```

```julia
julia> a = [0, 2, 1, 0, 0, 1];
```

```julia
julia> o=[1,3,4,2];
```

```julia

julia>  feynman_integral_branch_type_order(F,a,o) 
128*q[2]^4*q[3]^2*q[6]^2
```

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
julia>     substitute(feynman_integral_degree_sum(F,8))
10246144*q[1]^16 + 3294720*q[1]^14 + 886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
```

# Example of graph with loop.

![alt text](img/graph_loop.png)

```julia
julia> G=FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)])
graph([(1, 1), (1, 2), (2, 3), (3, 1)])
```

```julia
julia> F=FeynmanIntegral(G);
```

```julia
julia> O=[1,2,3]  
3-element Vector{Int64}:
 1
 2
 3
```

```julia
julia> a=[ 2,  0, 0, 1]
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
julia> substitute(feynman_integral_sum(F, 8))
20640*q[1]^16 + 9996*q[1]^14 + 4320*q[1]^12 + 1650*q[1]^10 + 456*q[1]^8 + 90*q[1]^6 + 6*q[1]^4
```


