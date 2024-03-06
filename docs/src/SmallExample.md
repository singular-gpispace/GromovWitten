```@meta
CurrentModule = GromovWitten
DocTestSetup = quote
  using GromovWitten
end
```

# First example

![alt text](img/caterpillar2.png)

To provide an example on how to use our package, we define a Feynman graph G from a list of edges:

```jldoctest mygraph
julia> G=FeynmanGraph([(1, 2), (1,2),(1,2)])
FeynmanGraph([(1, 2), (1, 2), (1, 2)])
```

We then define Feynman integral which includes polynomial  with all variables required by our implementation:

```jldoctest mygraph
julia> F=FeynmanIntegral(G)
FeynmanIntegral(FeynmanGraph([(1, 2), (1, 2), (1, 2)]), Dict{Symbol, Dict{Vector{Int64}, Nemo.QQMPolyRingElem}}(), (Multivariate polynomial ring in 7 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2]], Nemo.QQMPolyRingElem[q[1], q[2], q[3]], Nemo.QQMPolyRingElem[z[1], z[2]]))
```

We compute then the Feynman Integral of degree 3.

```jldoctest mygraph
julia> f = feynman_integral_degree(F, 3)
24*q[1]^6 + 20*q[1]^4*q[2]^2 + 20*q[1]^4*q[3]^2 + 20*q[1]^2*q[2]^4 + 20*q[1]^2*q[3]^4 + 24*q[2]^6 + 20*q[2]^4*q[3]^2 + 20*q[2]^2*q[3]^4 + 24*q[3]^6
```
