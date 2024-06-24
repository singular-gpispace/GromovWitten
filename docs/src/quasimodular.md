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
I(q)=\sum_{\substack{i,j,k \in \mathbb{N}_0 \\ 2i+4j+6k=6g-6}} b_{i,j,k} E_2^i E_4^j E_6^k
```

## Example

### Caterpillar genus 3

Consider the Caterpillar 3 graph

![alt text](img/Cartepillar3.png)

```jldoctest quasi
julia> G = FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
```

We then define the FeynmanIntegral type.

```jldoctest quasi
julia> F=FeynmanIntegral(G)
FeynmanIntegral(FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)]), Dict{Symbol, Dict{Vector{Int64}, Nemo.QQMPolyRingElem}}(), (Multivariate polynomial ring in 14 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3], x[4]], Nemo.QQMPolyRingElem[q[1], q[2], q[3], q[4], q[5], q[6]], Nemo.QQMPolyRingElem[z[1], z[2], z[3], z[4]]))
```

We compute the  sum of all Feynman Integral of degree up to 6.

```jldoctest quasi
julia> Iq=substitute(feynman_integral_degree_sum(F,6))
886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
```

To express the Feynman Integral Iq in term of Eisenstein series $E_2, E_4, E_6$, we define a polynomial ring in q

```jldoctest quasi
julia> R,q=polynomial_ring(QQ,["q"])
(Multivariate polynomial ring in 1 variable over QQ, Nemo.QQMPolyRingElem[q])
```

```jldoctest quasi
julia> Iqq=886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
886656*q^12 + 182272*q^10 + 25344*q^8 + 1792*q^6 + 32*q^4
```

we compute  quasimodular form of Iq :

```jldoctest quasi
julia> quasimodular_form(Iqq,12)
(1//93312, -3*E2^6 + 6*E2^4*E4 + 4*E2^3*E6 - 3*E2^2*E4^2 - 12*E2*E4*E6 + 4*E4^3 + 4*E6^2)
```

### Star graph $K_ {1,3}$

Consider the Star graph $K_ {1,3}$ .

![alt text](img/star_graph.png)

We want to find the quasimodular form of the sum of Feynman Integral  $Iq$ up to degree 6.

We define first a polynomial ring in one variable $q$.

```jldoctest form
julia> R,q=polynomial_ring(QQ,["q"])
(Multivariate polynomial ring in 1 variable over QQ, Nemo.QQMPolyRingElem[q])
```

The sum of Feynman Integral up to degree 6 is:

```jldoctest form
julia> Iq= 843264*q[1]^12 + 165888*q[1]^10 + 20736*q[1]^8 + 1152*q[1]^6 
843264*q^12 + 165888*q^10 + 20736*q^8 + 1152*q^6
```

The quasimodular form of $Iq$  is then:

```jldoctest form
julia> quasimodular_form(Iq,12)
(1//20736, -E2^6 + 3*E2^4*E4 - 3*E2^2*E4^2 + E4^3)
```
