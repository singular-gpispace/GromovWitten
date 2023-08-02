```@meta
CurrentModule = GromovWitten
DocTestSetup = quote
  using GromovWitten
end
```

# First example

![alt text](img/caterpillar2.png)

To provide an example on how to use our package, we define a graph G from a list of edges:

```jldoctest graph
julia> G=graph([(1, 2), (1,2),(1,2)])
graph([(1, 2), (1, 2), (1, 2)])
```

jldoctest graph
