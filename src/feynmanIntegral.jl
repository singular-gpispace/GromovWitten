function specificFeynmanIntegral(R::MPolyRing,x::Vector,q::Vector, G::graphe,a::Vector{Int64})
    #R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G))
    ee=Edge.(G.edge)

    N = sum(a)
    f = flip(G, a)
   p=0
   for i in 1:length(f)
      tmp=1
       for j in 1:length(f[i][2])
           if f[i][2][j] < 1
              if  f[i][2][j] == -1
                   tmp = tmp * constterm(x[src(ee[j])], x[dst(ee[j])], N)
               else
                   tmp = tmp * constterm(x[dst(ee[j])], x[src(ee[j])], N)
               end
           else
               tmp = tmp * proterm(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)
           end 
       end
      p=p+f[i][1]*coefterm(G,tmp,N)
           
   end
   return p
   
end


function feynmanIntegral(R::MPolyRing,x::Vector,q::Vector,G::graphe,d::Integer)
     ee=Edge.(G.edge)
     a=partition(length(ee),d) 
    sum=0
    for ai in a
        sum=sum+specificFeynmanIntegral(R,x,q,G,ai) 
    end
    return sum
end
function feynmanIntegralSum(R::MPolyRing,x::Vector,q::Vector,G::graphe,d::Integer)
    res=0
   for i in 1:d
        res+=feynmanIntegral(R,x,q,G,i)
    end
    return res
end

function subt(R::MPolyRing,x::Vector,q::Vector,p::fmpq_mpoly)

    if is_constant(R(p))
        return 0
    else
        r=gens(R)
        s=findfirst(isequal(q[1]),r)

       r[s:end] .= q[1]
    end
        return evaluate(p,r)
end