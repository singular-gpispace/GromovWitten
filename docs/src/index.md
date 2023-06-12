```@meta
CurrentModule = tropicalfeynman
```

# tropicalfeynman

Documentation for [tropicalfeynman](https://github.com/singular-gpispace/tropicalfeynman).

```@index

```

```@autodocs
Modules = [tropicalfeynman]
```

```
To run:
- pull the package in your local disk
In the local folder where the package is located, type in the terminal.
- julia --project
-  using tropicalfeynman
-  ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
-   G=graphe(ve)
-   R,x,q=polynomialring(G,"x","q")
-  a=[0,2,1,0,0,1]
-  specificFeynmanIntegral(R,x,q,G,a)
-  feynmanIntegral(R,x,q,G ,4)
-  subt(R,x,q,feynmanIntegral(R,x,q,G ,4))

```
