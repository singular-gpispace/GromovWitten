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

### Caterpillar genus 3

Consider the Caterpillar 3 graph

![alt text](img/Cartepillar3.png)

```jldoctest quasi
julia> G = graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
```

We compute the  sum of all Feynman Integral of degree up to 6.

```jldoctest quasi
julia> Iq=substitute(feynman_integral_degree_sum(G,6))
886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
```

To express the Feynman Integral Iq in term of Eisenstein series $E_2, E_4, E_6$.

We define a polynomial ring in $E2, E4, E6$

```jldoctest quasi
julia> S,E2,E4,E6=@polynomial_ring(QQ,E2,E4,E6)
(Multivariate polynomial ring in 3 variables over QQ, E2, E4, E6)
```

The quasimodular form of Iq is

```jldoctest quasi
julia> quasimodular_form(E2,E4,E6,q,Iq,12)
(1//93312, -3*E2^6 + 6*E2^4*E4 + 4*E2^3*E6 - 3*E2^2*E4^2 - 12*E2*E4*E6 + 4*E4^3 + 4*E6^2)
```

### Star graph $K_ {1,3}$

Consider the Star graph $K_ {1,3}$ .

![alt text](img/star_graph.png)

We want to find the quasimodular form of the sum of Feynman Integral  $Iq$ up to degree 6.

The sum of Feynman Integral up to degree 6 is:

```jldoctest form
julia> Iq= 843264*q^12 + 165888*q^10 + 20736*q^8 + 1152*q^6 
843264*q^12 + 165888*q^10 + 20736*q^8 + 1152*q^6
```

We define a polynomial ring in $E2, E4, E6$

```jldoctest form
julia> S,E2,E4,E6=@polynomial_ring(QQ,E2,E4,E6)
(Multivariate polynomial ring in 3 variables over QQ, E2, E4, E6)
```

The quasimodular form of $Iq$  is then:

```jldoctest form
julia> quasimodular_form(E2,E4,E6,q,Iq,12)
(1//20736, -E2^6 + 3*E2^4*E4 - 3*E2^2*E4^2 + E4^3)
```
