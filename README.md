# GromovWitten

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://singular-gpispace.github.io/GromovWitten/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://singular-gpispace.github.io/GromovWitten/

[ga-img]: https://github.com/singular-gpispace/GromovWitten/actions/workflows/CI.yml/badge.svg?branch=main
[ga-url]: https://github.com/singular-gpispace/GromovWitten/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/singular-gpispace/GromovWitten/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/singular-gpispace/GromovWitten

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ga-img]][ga-url] [![][codecov-img]][codecov-url] |

# GromovWitten

The package GromovWitten computes generating series for tropical Hurwitz numbers of elliptic curves via mirror symmetry and Feynman integrals, and thus, by a correspondence theorem, Hurwitz numbers in the sense algebraic geometry. Generalizations of the method also allow for the computation of Gromov-Witten invariants for ellptic curves, and are also implemented in the package. GromovWitten is based on the computeralgebra system OSCAR and is provided as a package for the Julia programming language.

# Installation

We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:

```bash
git clone https://github.com/singular-gpispace/GromovWitten.git
```

Navigate to  GromovWitten  folder execute the following command:

```bash
julia --project
```

This will activate the environment for our package. In Julia install missing packages:

```bash
import Pkg; Pkg.instantiate()
```

and load our package. On the first run this may take some time.

```bash
using GromovWitten  
```


# Example of graph without vertex contribution and loop.

![alt text](docs/src/img/Cartepillar3.png)

```julia
julia> G = feynman_graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
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


# Caching Feynman Integral

We can try to catch the previous result in a table. 
To do that, we define `feynman_integral_branch_type_cache` , `feynman_integral_degree_cache` and  `feynman_integral_degree__sum_cache`. The last one returns the univariable polynomial of Feynman integral sum up to degree $d$.
We define first the Feynman graph, 
```julia
julia> G = feynman_graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
```

```julia
julia> F=FeynmanIntegral(G)
```
we cache the Feynman Integral.
```julia
julia> feynman_integral_degree_cache(F, 3);
```
To diplay the previous caching 
```julia
julia> F.integral_cache[:degree]
Dict{Vector{Int64}, QQMPolyRingElem} with 4 entries:
  [8]  => 906376*q[1]^16 + 76832*q[1]^14*q[2]^2 + 76832*q[1]^14*q[3]^2 + 76832*…
  [3]  => 288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 +…
  [5]  => 20000*q[1]^10 + 2592*q[1]^8*q[2]^2 + 2592*q[1]^8*q[3]^2 + 2592*q[1]^8…
  [10] => 5465008*q[1]^20 + 350352*q[1]^18*q[2]^2 + 350352*q[1]^18*q[3]^2 + 350…
```
We get the following table of comparison.
|Type|degree|first_run|second_run|
|:----|:----|:----|:----|
|degree|[3] => 288*q[1]^6…|0.044741 s(17.840 MiB)|0.000021 s ( 1.609 KiB)|
|degree|[5] => 20000*q[1]^10…| 0.209154 s( 129.200 MiB)|0.000021 s ( 1.609 KiB)|
|degree| [8] => 906376*q[1]^16 …|1.768598 s( 1.113 GiB)|1.000021 s ( 1.609 KiB)|
|degree|[10] => 5465008*q[1]^20 +…|5.453744 s(3.425 GiB)|0.000033 s ( 1.609 KiB)|


simillary for The univariable polynomial, we have 
```julia
julia> feynman_integral_degree_sum_cache(F, 3)
1792*q[1]^6 + 32*q[1]^4
```
To diplay the previous caching 
```julia
julia> F.integral_cache[:sum]
Dict{Vector{Int64}, QQMPolyRingElem} with 5 entries:
  [8]  => 10246144*q[1]^16 + 3294720*q[1]^14 + 886656*q[1]^12 + 182272*q[1]^10 …
  [3]  => 1792*q[1]^6 + 32*q[1]^4
  [5]  => 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
  [11] => 145337600*q[1]^22 + 66497472*q[1]^20 + 27353088*q[1]^18 + 10246144*q[…
  [10] => 66497472*q[1]^20 + 27353088*q[1]^18 + 10246144*q[1]^16 + 3294720*q[1]…
```

|Type|degree|normal_time|cache_time|
|:----|:----|:----|:----|
|sum| [3] => 1792*q[1]^6…|  0.028567 s ( 23.850 MiB) |  0.029917 s (23.852 MiB)sum|[5] => 182272*q[1]^10…|  0.337313 s (203.303 MiB) |  0.293910 s (179.465 MiB)sum|[8] => 10246144*q[1]^16  …|  3.529455 s  (2.169 GiB)| 3.210641 s (1.970 GiB|
|sum|[10] => 66497472*q[1]^20  +…| 12.292280 s (7.607 GiB| 8.813801 s ( 5.439 GiB)|
|sum |[11] => 145337600*q[1]^22+…|21.427653 s (13.190 )|9.026728 s (5.583 GiB)|
