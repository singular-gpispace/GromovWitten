```@meta
CurrentModule = GromovWitten
DocTestSetup = quote
using GromovWitten
end
```

# Quasimodular

In this module, we compute the solution of the system $Ax=b$, where $A$ is a matrix from homogeneous Eisenstein series $E_2, E_4, E_6$ and $b$
from the Feynman Integral $I(q)$
The solution is of the form (factor, coefficients) where coefficients is a vector of rationals numbers.
Given a Feynman Integral

```math
I(q)=\sum_{n=1}^{d} a_i q^{d},
```

we compute the coefficients $b_{i,j,k}$ such that

```math
I(q)=\sum_{\substack{i,j,k \in \mathbb{N}_0 \\ 2i+4j+6k=d}} b_{i,j,k} E_2^i E_4^j E_6^k
```

## Example

### Graph with vertex contribution

Consider the Graph with vertex contribution

![alt text](img/graph_with_vertex1.png)

```jldoctest quasi
julia> G = FeynmanGraph( [(1, 2), (2, 3), (3, 1)])
FeynmanGraph([(1, 2), (2, 3), (3, 1)])
```

We then define the FeynmanIntegral type.

```jldoctest quasi
julia> F=FeynmanIntegral(G)
FeynmanIntegral(FeynmanGraph([(1, 2), (2, 3), (3, 1)]), Dict{Symbol, Dict{Vector{Int64}, Nemo.QQMPolyRingElem}}(), (Multivariate polynomial ring in 9 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3]], Nemo.QQMPolyRingElem[z[1], z[2], z[3]]))
```

We compute the  sum of all Feynman Integral of degree up to $m=\text{number\_of\_monomials}(\text{weightmax}) $
where weightmax$=2(r+\sum_{i=1}^{n} g_i)$ with $r$ the number of edges and $g_i$ satifying  $h^1(\Gamma)+\sum_{i=1}^{n} g_i=g$ .
Suppose $gg=[1,0,0]$, we have r=3; so
weightmax$=2(3+1)=8$

```jldoctest quasi
julia> weightmax=8;
```

```jldoctest quasi
julia> m = number_of_monomial(weightmax)
10  

```


We computed then the Feynman Integral sum up to the degree $m=10$

```jldoctest quasi
julia> Iq=substitute(feynman_integral_degree_sum(F, m,g))
56250*q[1]^20 + 121581//4*q[1]^18 + 18480*q[1]^16 + 8330*q[1]^14 + 4428*q[1]^12 + 3075//2*q[1]^10 + 556*q[1]^8 + 117*q[1]^6 + 15*q[1]^4 + 1//4*q[1]^2
```

We can now express the Feynman Integral Iq in term of Eisenstein series $E_2, E_4, E_6$ by call the

we compute  quasimodular form of Iq :

```jldoctest quasi
julia> quasimodularity_form(Iq,weightmax)
(1//6912, E2^3 + 2*E2^2*E4 - 3*E2*E4 - 4*E2*E6 + 2*E4^2 + 2*E6)
```

### Graph with loop at the vertex 1

Consider the Graph with loop at the vertex 1 .

![alt text](img/graph_loop.png)

```jldoctest quasi
julia> G = FeynmanGraph( [(1, 1),(1, 2), (2, 3), (3, 1)])
FeynmanGraph( [(1, 2), (2, 3), (3, 1)])
```

We then define the FeynmanIntegral type.

```jldoctest quasi
julia> F=FeynmanIntegral(G)
FeynmanIntegral(FeynmanGraph( [(1, 1),(1, 2), (2, 3), (3, 1)]), Dict{Symbol, Dict{Vector{Int64}, Nemo.QQMPolyRingElem}}(), (Multivariate polynomial ring in 10 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3], q[4]], Nemo.QQMPolyRingElem[z[1], z[2],z[3]]))
```

We compute the  sum of all Feynman Integral of degree up to $m=\text{number\_of\_monomials}(\text{weightmax}) $.
Here $gg=[0,0,0]$, and r=4; so
weightmax$=2(4+0)=8$.

```jldoctest quasi
julia> weightmax=8;
julia> m=number_of_monomial(weightmax)
10
```

We computed then the Feynman Integral sum up to the degree $m=10$

```jldoctest quasi
julia> Iq=substitute(feynman_integral_degree_sum(F, m))
67500*q[1]^20 + 36774*q[1]^18 + 20640*q[1]^16 + 9996*q[1]^14 + 4320*q[1]^12 + 1650*q[1]^10 + 456*q[1]^8 + 90*q[1]^6 + 6*q[1]^4
```

We can now express the Feynman Integral Iq in term of Eisenstein series $E_2, E_4, E_6$ by call the

we compute  quasimodular form of Iq :

```jldoctest quasi
julia> quasimodularity_form(Iq,weightmax)
(1//6912, E2^4 - E2^3 - 3*E2^2*E4 + 3*E2*E4 + 2*E2*E6 - 2*E6)
```

