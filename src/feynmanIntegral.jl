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
function specificFeynmanIntegralV(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector, G::graphe, a::Vector{Int64};aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipV(G, a)
    p=0
     sz=1
    for k in 1:length(l)
        sz=sz*InvSfunction(z[k],aa)
    end
    for i in 1:length(f)
        tmp=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * consttermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N)

                elseif f[i][2][j] == 0
                        tmp = tmp * consttermV(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N)

                elseif f[i][2][j] == -2
                        tmp = tmp *looptermV(z[j],q[j],aa,a[j])

                else 
                        tmp = tmp * protermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N)
                end 
        end
        p=p+f[i][1]*coeftermX(G,tmp,l,N)
        
        if i==length(f)
            pr=sz*p
            p=coefterm2Z(G,pr,g)
        end

            
    end
    return p
end

function feynmanIntegralV(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector, G::graphe,d::Integer; aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegralV(R,x,q,z,G,ai;aa,l,g) 
    end
    return sum
end 
function feynmanIntegralSumV(R::Nemo.FmpqMPolyRing,x::Vector,q::Vector,z::Vector,G::graphe,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    res=0
    for i in 1:d
        res+=feynmanIntegralV(R,x,q,z,G,i;aa,l,g)
    end
    return res
end
   
function subtV(p::fmpq_mpoly)
coeffs_dict = coeffs(p)
coeffs_array = collect(values(coeffs_dict))
return sum(coeffs_array)
end