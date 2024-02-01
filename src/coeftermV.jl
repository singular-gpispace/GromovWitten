###############################################################################
#                                                                             #
#    coeftermV.jl : compute differents  function for vertex contributions.    #
#                                                                             #
###############################################################################
@doc raw"""
    sfunction(z::QQMPolyRingElem,k::Int64)

**Note**:The function sfunction(z,k) takes account vertex contributions. 
```math
S(z, {aa}) = \sum_{n = 0}^{{aa}} \dfrac{2^{- 1 - n} (1 + (- 1)^n)
}{(n + 1) !} z^n = \sum_{n = 0}^{{aa}} \dfrac{{2^{- 2 n}} }{(2 n + 1) !}
z^n, {aa} \rightarrow \infty
```
# Examples

```julia
julia> R,x=polynomial_ring(QQ,:x=>1:1); # using Nemo
julia> sfunction(x[1],4)

1//92897280*z[1]^8 + 1//322560*z[1]^6 + 1//1920*z[1]^4 + 1//24*z[1]^2 
```
"""
function sfunction(x::QQMPolyRingElem,k::Int64)
    su=0
   for n in 0:k
       su=su+(1//BigInt(2^(2n)))//(factorial(BigInt((2*n+1))))*(x)^(2*n)
       
   end
   return su
end

@doc raw"""
    inv_sfunction(z::QQMPolyRingElem,aa::Int64)

returns the inverse sfunction
```math
 \frac{1}{S(z,aa)}=\frac{z}{2 Sinh(z/2)}= 
 \sum_{n =  0}^{aa} \left( \left\{\begin{array}{ll}
     1 & \text{if} && n = 1\\
     - \frac{- 2^n (- 2 + 2^n)}{n!} B_n & (n > = 1 && (- 1
     + n)\mod 2 = 1)
   \end{array}\right. \right) z^n 
```
   Where $B_n$ is Bernoulli number and ${aa} \rightarrow \infty$.
# Examples
```julia
julia> R,x=polynomial_ring(QQ,:x=>1:1); # using Nemo
julia> inv_sfunction(x[1],4)
7//5760*x[1]^4 - 1//24*x[1]^2 + 1
```
"""
function inv_sfunction(z::QQMPolyRingElem,k::Int64)
    su=0
    for n in 0:k+1
        su=su-(((1//2^(n))*(-2 + 2^n)*bernoulli(n))//(factorial(n)))*z^n
        
    end
    return su
end
@doc raw"""
    loopterm( q::QQMPolyRingElem, a::Integer)

returns loop contribution with zero genus gi at a vertex i. 
"""

function loopterm( q::QQMPolyRingElem, a::Integer)
    p=0
   if a==0
       return p
   else 
       for w in 1:a
           if a%w==0
               p = p + w*q^(2*a)
           end
       end
   end
       return p
end
@doc raw"""
    loopterm( z::QQMPolyRingElem, q::QQMPolyRingElem, aa::Integer, a::Integer)

returns loop contribution with nonzero genus gi at a vertex i. 
"""
function loopterm( z::QQMPolyRingElem,q::QQMPolyRingElem, aa::Integer, a::Integer)
    p=0
   if a==0
       return p
   else 
       for w in 1:a
           if a%w==0
               S1=sfunction(w*z,aa)
               p = p + S1*S1*w*q^(2*a)
           end
       end
   end
       return p
end
#=function coefterm2Z( x::Vector, q::Vector,z::Vector,hp::QQMPolyRingElem,g::Vector) 
    g=2 .* g
    if hp==0
        return 0
    else
        hp=coeff(hp,z,g)
    end
    return hp
end
function coeftermX( x::Vector, q::Vector,z::Vector,G::graph ,p::QQMPolyRingElem,d::Integer;l=zeros(Int,nv(G)))
    ee=edges(G)
    G=DiGraph(edges(G))
    L=zeros(Int,nv(G))
    for ev in ee
        if src(ev) == dst(ev)

        else
            L[src(ev)]=L[src(ev)]+d
            L[dst(ev)]=L[dst(ev)]+d
        end
    end
    L=L .+l
    p=coeff(p,x,L)
    return p
end =#
function lis(G::FeynmanGraph,d::Int64,l::Vector{Int64})
    ee = edges(G)
    #G=DiGraph(edges(G))
    L=zeros(Int,nv(G))
     @inbounds for ev in ee
        if src(ev) != dst(ev)
            L[src(ev)] += d
            L[dst(ev)] += d
        end
    end
               L=L .+l
   return L
end
@doc raw"""
     filter_term(p::Union{QQMPolyRingElem, Int64}, variables::Vector{QQMPolyRingElem}, s::Vector{Int64})
 
 replaces all terms of the polynomial p with zero whenever the variables raised to a power of s1 exceed the specified power s.
  
# Examples
```julia  
julia> p=8*q[1]^6*q[2]^2 + 8*q[1]^6*q[3]^2 + 8*q[1]^6*q[4]^2 + 54*q[1]^4*q[2]^4 + 18*q[1]^4*q[2]^2*q[3]^2 + 18*q[1]^4*q[2]^2*q[4]^2 
+ 54*q[1]^4*q[3]^4 + 18*q[1]^4*q[3]^2*q[4]^2 + 54*q[1]^4*q[4]^4 + 56*q[1]^2*q[2]^6 + 6*q[1]^2*q[2]^4*q[3]^2 + 6*q[1]^2*q[2]^4*q[4]^2 
+ 6*q[1]^2*q[2]^2*q[3]^4 + 12*q[1]^2*q[2]^2*q[3]^2*q[4]^2 + 6*q[1]^2*q[2]^2*q[4]^4 + 56*q[1]^2*q[3]^6 + 6*q[1]^2*q[3]^4*q[4]^2 + 6*q[1]^2*q[3]^2*q[4]^4 + 56*q[1]^2*q[4]^6
```
we replace all term in $p$  with `q[1]^a*q[2]^b*q[3]^c > q[1]*q[2]*q[3]` by zero,this means all power $(a,b,c)>(2,2,2)$

```julia  
julia> filter_term(p,[q[1],q[2],q[3]],[2,2,2])
12*q[1]^2*q[2]^2*q[3]^2*q[4]^2 + 6*q[1]^2*q[2]^2*q[4]^4 + 6*q[1]^2*q[3]^2*q[4]^4 + 56*q[1]^2*q[4]^6
```
also
```julia  
julia>  filter_term(p,[q[1],q[2]],[2,2])
6*q[1]^2*q[2]^2*q[3]^4 + 12*q[1]^2*q[2]^2*q[3]^2*q[4]^2 + 6*q[1]^2*q[2]^2*q[4]^4 + 56*q[1]^2*q[3]^6 + 6*q[1]^2*q[3]^4*q[4]^2 + 6*q[1]^2*q[3]^2*q[4]^4 + 56*q[1]^2*q[4]^6 + q[1]

julia> filter_term(p,q[1],1)
q[1]
```
 """
 function filter_term(pols::Union{QQMPolyRingElem,QQPolyRingElem, Int64}, variables::Union{Vector{QQPolyRingElem},Vector{QQMPolyRingElem}, QQMPolyRingElem,QQPolyRingElem}, power::Union{Vector{Int64}, Int64})
    if typeof(variables) <: QQMPolyRingElem
        T = parent(variables)
        pols = T(pols)
        gensR = gens(T)
        position = findfirst(var -> var == variables, gensR)
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
    else
        T = parent(variables[1])
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
end

