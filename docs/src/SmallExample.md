```@meta
CurrentModule = GromovWitten
```

# Examples Hurwitz numbers.

![alt text](img/caterpillar2.png)


To provide an example on how to use our package, we define a graph G from a list of edges:

```julia
julia> v=[(1, 2), (1,2),(1,2)]
```

```julia
julia> G=graph([(1, 2), (1,2),(1,2)])
```

We then define a polynomial ring with all variables required by our implementation:

```julia
julia> R, x, q = polynomial_ring(G, "x", "q")
(Multivariate polynomial ring in 10 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3],x[4]], QQMPolyRingElem[q[1], q[2],q[3], q[4], q[5], q[6]])
```
