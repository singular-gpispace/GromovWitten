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
julia> R,x=@polynomial_ring(QQ,x[1:1]); # using Nemo
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
julia> R,x=@polynomial_ring(QQ,x[1:1]); # using Nemo
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
               p = p + w*q^(a)
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
               p = p + S1*S1*w*q^(a)
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
    ee=Edge.(G.edge)
    G=DiGraph(Edge.(G.edge))
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
function lis(G::graph,d::Int64,l::Vector{Int64})
    ee = Edge.(G.edge)
    #G=DiGraph(Edge.(G.edge))
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
julia> p=feynman_integral_degree(x,q,G,4)
 8*q[1]^3*q[2] + 8*q[1]^3*q[3] + 8*q[1]^3*q[4] + 54*q[1]^2*q[2]^2 + 18*q[1]^2*q[2]*q[3] 
 + 18*q[1]^2*q[2]*q[4] + 54*q[1]^2*q[3]^2 + 18*q[1]^2*q[3]*q[4] + 54*q[1]^2*q[4]^2 
 + 56*q[1]*q[2]^3 + 6*q[1]*q[2]^2*q[3] + 6*q[1]*q[2]^2*q[4] + 6*q[1]*q[2]*q[3]^2 
 + 12*q[1]*q[2]*q[3]*q[4] + 6*q[1]*q[2]*q[4]^2 + 56*q[1]*q[3]^3 + 6*q[1]*q[3]^2*q[4] 
 + 6*q[1]*q[3]*q[4]^2 + 56*q[1]*q[4]^3
```

we replace all term in $p$  with `q[1]^a*q[2]^b*q[3]^c > q[1]*q[2]*q[3]` by zero,this means all power $(a,b,c)>(1,1,1)$

```julia  
julia> filter_term(p,[q[1],q[2],q[3]],[1,1,1])
 12*q[1]*q[2]*q[3]*q[4] + 6*q[1]*q[2]*q[4]^2 + 6*q[1]*q[3]*q[4]^2 + 56*q[1]*q[4]^3
```
 """
function filter_term(p::Union{QQMPolyRingElem, Int64}, variables::Vector{QQMPolyRingElem}, s::Vector{Int64})
    T = parent(variables[1])
  
    gensR = gens(T)
    position = [findfirst(var -> var == vi, gensR) for vi in variables]
    result = zero(p)
    
    d = Vector{Int}(undef, length(variables))  # Preallocate d
    @inbounds for term in terms(T(p))
        for j in 1:length(variables)
            po = position[j]
            d[j] = degree_fmpz(term, po)
        end
        if all(d .<= s)
            result += term
        end
    end
    
    return result
end
