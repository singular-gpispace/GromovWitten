function express_as_eisenstein_series(E2,E4,E6,n::Int)
    expressions = []
    for e2 in 0:n
        for e4 in 0:n
            for e6 in 0:n
                d = 2*e2 + 4*e4 + 6*e6
                if d == n && (e2 > 0 || e4 > 0 || e6 > 0)
                    term = 1
                    if e2 != 0
                        term *= E2^e2
                    end
                    if e4 != 0
                        term *= E4^e4
                    end
                    if e6 != 0
                        term *= E6^e6
                    end
                    push!(expressions, term)
                end
            end
        end
    end
    return expressions
end

#=function filter_term(pols::Union{QQMPolyRingElem, Int64}, variables::Vector{QQMPolyRingElem}, power::Vector{Int64})
    T = parent(variables[1])
    if typeof(pols) <: Integer
        pols = T(pols)
    end
    gensR = gens(T)
    position = [findfirst(var -> var == vi, gensR) for vi in variables]
    result = zero(pols)
    d = Vector{Int}(undef, length(variables))  # Preallocate d
    @inbounds for term in terms(T(pols))
        for j in 1:length(variables)
            po = position[j]
            d[j] = degree_fmpz(term, po)
        end
        if all(d .<= power)
            result += term
        end
    end
    
    return result
end
=#
@doc raw"""
    filter_vector(polyvector::Vector{QQMPolyRingElem}, variables::Union{Vector{QQMPolyRingElem}, QQMPolyRingElem}, power::Union{Vector{Int64}, Int64})

returns the number filter_term of function in a vector of polynomial.
# Example
```julia
julia> R,q=@polynomial_ring(QQ,q)

julia> v=[-122976*q^6 - 16632*q^4 - 504*q^2 + 1, -645120*q^12 - 691200*q^10 - 339840*q^8 - 62496*q^6 - 3672*q^4 + 216*q^2 + 1, -884736*q^18 - 1990656*q^16 - 2156544*q^14 - 1340928*q^12 - 497664*q^10 - 95040*q^8 - 3744*q^6 + 1512*q^4 - 72*q^2 + 1]

julia> filter_vector(v,q,6)
3-element Vector{QQMPolyRingElem}:
 -122976*q^6 - 16632*q^4 - 504*q^2 + 1
 -62496*q^6 - 3672*q^4 + 216*q^2 + 1
 -3744*q^6 + 1512*q^4 - 72*q^2 + 1
```
"""
function filter_vector(polyvector::Vector{QQMPolyRingElem}, variables::Union{Vector{QQMPolyRingElem}, QQMPolyRingElem}, power::Union{Vector{Int64}, Int64})
    result = Vector{QQMPolyRingElem}()
    for pols in polyvector
        push!(result, filter_term(pols, variables, power))
    end
    return result
end
@doc raw"""
    sum_of_divisor_powers(n::Int, k::Int)

$σ_k(n)$  returns the sum of the $k^{th}$ powers of divisors of $n$.
        it returns also the number of divisors $d(n)$  of $n$ for $k=0$. 
```julia
julia>sum_of_divisor_powers(6,0) # number of divisors of 6 $d(6)$
4

sum of divisors of 6.
julia> sum_of_divisor_powers(6,1)
12
```
"""
function sum_of_divisor_powers(n::Int, k::Int)
    if n==0
        error(" First Argument must be non-zero")
    end
    divs = divisors(n)
    div_sum = sum(div^k for div in divs)
    return div_sum
end
# Compute the Eisenstein series E2(q[1]) up to a specified number of terms

@doc raw"""
     eisenstein_series(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}}, order::Int,k::Int)                  

Return the expansion of the  weight  Eisenstein series k with fixed order.
For a fixed order $m$, we compute $$E_k = 1 - \frac{2 k}{ B_k}  \sum_{d=1}^{m} \sigma_{k-1}(d)  q^{2 d}$$
```julia
julia> eisenstein_series(q,5,2)  #E2
-144*q^10 - 168*q^8 - 96*q^6 - 72*q^4 - 24*q^2 + 1

julia> eisenstein_series(q,5,4) #E4
30240*q^10 + 17520*q^8 + 6720*q^6 + 2160*q^4 + 240*q^2 + 1
```
"""
function eisenstein_series(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}}, num_terms::Int,k::Int)                  
    if k % 2 != 0
        error("input k must be even in eisenstein_series(q, num_terms,k)")
    end
    if typeof(q) == Vector{QQMPolyRingElem}
            return  e2 = 1 - ((((2*k)// bernoulli(k))) )* sum(sum_of_divisor_powers(d, k-1) * q[1]^(2*d) for d in 1:num_terms)
    else
            return  e2 = 1 - ((((2*k)// bernoulli(k))) )* sum(sum_of_divisor_powers(d, k-1) * q^(2*d) for d in 1:num_terms)
    end
end

function express_as_powers(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}}, max_degree::Int)
    if max_degree % 2 != 0
        error("input must be even")
    end
    result = Vector{QQMPolyRingElem}()
    nb=Int64(max_degree//2)
    E6=eisenstein_series(q, nb,6)
    E4=eisenstein_series(q, nb,4)
    E2=eisenstein_series(q, nb,2)
    for e2 in 0:max_degree
        for e4 in 0:max_degree
            for e6 in 0:max_degree
                degree = 2*e2 + 4*e4 + 6*e6
                if degree == max_degree
                    push!(result,  (E2^e2) * (E4^e4) * (E6^e6))
                end
            end
        end
    end
    return result
end
@doc raw"""
     polynomial_to_matrix(vect::Vector{QQMPolyRingElem})

returns a matrix from a given vector of polynomials with same degree. The returned matrix is of type QQMatrix.

```julia
julia> vp=[ -122976*q^6 - 16632*q^4 - 504*q^2 + 1, -62496*q^6 - 3672*q^4 + 216*q^2 + 1,-3744*q^6 + 1512*q^4 - 72*q^2 + 1]

julia> polynomial_to_matrix(vp)
[      1        1       1]
[      0        0       0]
[   -504      216     -72]
[      0        0       0]
[ -16632    -3672    1512]
[      0        0       0]
[-122976   -62496   -3744]
```
"""
function polynomial_to_matrix(vect::Vector{QQMPolyRingElem})
    A = matrix(QQ, hcat([collect(coefficients_of_univariate(poly)) for poly in vect]...))
    return A
end
@doc raw"""
     polynomial_to_matrix(vect::Vector{QQMPolyRingElem})

returns a matrix from a given   polynomial. The returned matrix is of type QQMatrix.

```julia
julia> Iq=25344*q^8 + 1792*q^6 +32q^4

julia> matrix_of_integral(Iq)
[    0]
[    0]
[    0]
[    0]
[   32]
[    0]
[ 1792]
[    0]
[25344]
```
"""
function matrix_of_integral(Iq::QQMPolyRingElem) 
    # Obtain the coefficients of Iq
    Q = coefficients_of_univariate(Iq)
        
    # Convert Q to a matrix representation
    Q_matrix = Matrix(transpose(reduce(hcat, Q)))

    # Convert the matrix to the QQ type so we can call can_with_solution.
    Q = matrix(QQ, Q_matrix)
        return Q
end

function solve_polynomial_system(A::QQMatrix, Q::QQMatrix)
    m, n = size(A)
    
    # Solve the polynomial system
    can_solve, x = can_solve_with_solution(A, Q)
    
    # Compute the common factor and scaled vector
    if can_solve
        common_factor = gcd(x...)
        coeff_vector = [val // common_factor for val in x]
        return  common_factor, coeff_vector
    else
        return "The system has no solution"
    end
end
#=@doc raw"""
     quasi_matrix(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}},Iq::QQMPolyRingElem, max_degree::Int64)

returns solution of the system $Ax=b$, where A is a matrix from homogeneous Eisenstein series $E_2, E_4, E_6$ and $b$
from the Feynman Integral $I(q)$
The solution is of the form (factor, coefficients) where coefficients is a vector of rationals numbers.
Given a Feynman Integral $$I(q)=\sum_{n=1}^{d} a_i q^{d}$$, we compute the coefficients $b_{i,j,k}$ such that 
$$I(q)=\sum_{i,j,k} b_{i,j,k} E_2^i E_4^j E_6^k$$

```julia
julia> Iq=886656*q^12 + 182272*q^10 + 25344*q^8 + 1792*q^6 +32q^4

julia> quasi_matrix(q,Iq,12)
(1//93312, QQFieldElem[4; 4; -12; -3; 4; 6; -3])
```
"""=#
function quasi_matrix(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}},Iq::QQMPolyRingElem, max_degree::Int64)
    if  total_degree(Iq)!=max_degree
        throw(DimensionMismatch("Degree of the polynomial must be equal to max_degree  $max_degree"))
    else  
        Evector=filter_vector(express_as_powers( q,max_degree),q,max_degree) 
        A=polynomial_to_matrix(Evector)
        Q=matrix_of_integral(Iq)
        return solve_polynomial_system(A,Q)
    end
end
@doc raw"""
     quasimodular_form(q, Iq, max_degree)

express the Feynman Integral polynomial $I(q)$ in terms of a polynomial in  $E_2, E_4, E_6$ 
This leads to 
$$I(q)=\sum_{i,j,k} b_{i,j,k} E_2^i E_4^j E_6^k$$

```julia
julia> R,q=@polynomial_ring(QQ,q)

julia> Iq=886656*q^12 + 182272*q^10 + 25344*q^8 + 1792*q^6 +32q^4
```
We define the polynomial ring in E2,E4,E6.

```julia
julia> S,E2,E4,E6=@polynomial_ring(QQ,E2,E4,E6)

julia> quasimodular_form(q,Iq,12)
(1//93312, -3*E2^6 + 6*E2^4*E4 + 4*E2^3*E6 - 3*E2^2*E4^2 - 12*E2*E4*E6 + 4*E4^3 + 4*E6^2)
```
"""
function quasimodular_form(E2,E4,E6,q, Iq, max_degree)
    fac, coef = quasi_matrix(q, Iq, max_degree)
    comb_result = express_as_eisenstein_series(E2,E4,E6,12)
    p=0
    for (i, term) in enumerate(comb_result)
        if coef[i] == 0
            continue
        end
        tmp=coef[i]*term
        p=p+tmp
    end
    p=fac,p
    
    return p
end

