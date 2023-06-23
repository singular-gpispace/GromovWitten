###############################################################################
#                                                                             #
#    coeftermV.jl : compute differents  function for vertex contributions.    #
#                                                                             #
###############################################################################
@doc raw"""
Sfunction(z::QQMPolyRingElem,k::Int64)
**Note**:The function Sfunction(z,k) takes account vertex contributions. 
$S (z, {aa}) = \sum_{n = 0}^{{aa}} \dfrac{2^{- 1 - n} (1 + (- 1)^n)
}{(n + 1) !} z^n = \sum_{n = 0}^{{aa}} \dfrac{{2^{- 2 n}} }{(2 n + 1) !}
z^n, {aa} \rightarrow \infty$
# Examples
```jldoctest
julia> Sfunction(z[1],4)
1//92897280*z[1]^8 + 1//322560*z[1]^6 + 1//1920*z[1]^4 + 1//24*z[1]^2 
````
"""
function Sfunction(z::QQMPolyRingElem,k::Int64)
    su=0
   for n in 0:k
       su=su+(1//(2^(2n)))//(factorial(2*n+1))*(z)^(2*n)
   end
   return su
end
function InvSfunction(z::QQMPolyRingElem,k::Int64)
    su=0
    for n in 0:k+1
        su=su-(((1//2^(n))*(-2 + 2^n)*bernoulli(n))//(factorial(n)))*z^n
        
    end
    return su
end
function looptermV( q::fmpq_mpoly, a::Integer)
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
function looptermV( z::fmpq_mpoly,q::fmpq_mpoly, aa::Integer, a::Integer)
    p=0
   if a==0
       return p
   else 
       for w in 1:a
           if a%w==0
               S1=Sfunction(w*z,aa)
               p = p + S1*S1*w*q^(a)
           end
       end
   end
       return p
end
function coefterm2Z(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector,hp::fmpq_mpoly,g::Vector) 
    g=2 .* g
    if hp==0
        return 0
    else
        hp=coeff(hp,z,g)
    end
    return hp
end
function coeftermX(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector,G::graph ,p::fmpq_mpoly,d::Integer;l=zeros(Int,nv(G)))
    r=gens(R)
    m=r[1:nv(G)]
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
    p=coeff(p,m,L)
    return p
end
