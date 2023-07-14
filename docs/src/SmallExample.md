```@meta
CurrentModule = GromovWitten
DocTestSetup = quote
  using GromovWitten
end
```

# First example

![alt text](img/caterpillar2.png)


To provide an example on how to use our package, we define a graph G from a list of edges:

```jldoctest 
julia> G=graph([(1, 2), (1,2),(1,2)])
 graph([(1, 2), (1, 2), (1, 2)])
```
We then define a polynomial ring with all variables required by our implementation:

```jldoctest
julia> R, x, q = polynomial_ring(G, "x", "q")
 (Multivariate polynomial ring in 5 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2]], Nemo.QQMPolyRingElem[q[1], q[2], q[3]])
```
```jldoctest
julia> f = feynman_integral_degree(x, q, G, 3)
24*q[1]^6 + 20*q[1]^4*q[2]^2 + 20*q[1]^4*q[3]^2 + 20*q[1]^2*q[2]^4 + 20*q[1]^2*q[3]^4 + 24*q[2]^6 + 20*q[2]^4*q[3]^2 + 20*q[2]^2*q[3]^4 + 24*q[3]^6
```
