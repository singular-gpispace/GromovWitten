```@meta
CurrentModule = GromovWitten
```

# GromovWitten

Documentation for [tropicalfeynman](https://github.com/singular-gpispace/tropicalfeynman).

```@index

```

```@autodocs
Modules = [GromovWitten]
```

```
To run:
- pull the package in your local disk
In the local folder where the package is located, type in the terminal.
- julia --project
-  using GromovWitten
-  ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
-   G=graph(ve)
-   R,x,q=polynomial_ring(G,"x","q")
-  a=[0,2,1,0,0,1]
-  feynman_integral_branchtype(x,q,G,a)
-  feynman_integral_degree(x,q,G ,4)
-  subt(R,x,q,feynman_integral_degree(x,q,G ,4))

```
