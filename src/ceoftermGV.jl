function Sfunction(x,k)
    su=0
   for n in 0:k
       su=su+(1//(2^(2n)))//(factorial(2*n+1))*(x)^(2*n)
       
   end
   return su
end
function InvSfunction(x,k)
    su=0
    for n in 0:k+1
        su=su-(((1//2^(n))*(-2 + 2^n)*Oscar.bernoulli(n))//(factorial(n)))*x^n
        
    end
    return su
end
function loopterm( z::fmpq_mpoly,q::fmpq_mpoly, aa::Integer, a::Integer)
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
function coeftermQ(G::graphe,hp::fmpq_mpoly,a::Vector) 
    nn=nv(G)+1
    n=ne(G)
    r=gens(R)
    m=length(r)
    if hp==0
        return 0
    else
        for i in nn:nn+n-1
            f=hp
            u=i-nn+1
            if a[i-nn+1]==0
            else
                m=coefficients(f,i)
                gg=a[i-nn+1]
                hp=m[gg+1] 
            end
        end
    end
    return hp
end
function coefterm2Z(G::graphe,hp::fmpq_mpoly,g::Vector) 
    nn=nv(G)+ne(G)+1
    n=nv(G)+ne(G)+nv(G)
    if hp==0
        return 0
    else
        for i in nn:n
            f=hp
            m=coefficients(f,i)
            gg=g[i-nn+1]
            hp=m[2*gg+1]
        end
    end
        return hp
end
function coeftermX(G::graphe ,p::fmpq_mpoly,l::Vector,d::Integer)
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
 for (i,li) in enumerate(L) 
     f=p
     m=coefficients(f,i)
     if size(m,1)<=li
         return 0 
     end
     p=m[li+l[i]+1,1]
 end
 return p
end
