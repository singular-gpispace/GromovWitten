function specificFeynmanIntegral(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector, G::graphe ,a::Vector{Int64} ;l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipV(G, a)
    p=0 
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * constterm(x[src(ee[j])], x[dst(ee[j])],N)

                elseif f[i][2][j] == 0
                        tmp = tmp * constterm(x[dst(ee[j])], x[src(ee[j])], N)
                else 
                        tmp = tmp * proterm(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)
                end 
        end
            #p=p+f[i][1]*tm
        p=p+f[i][1]*coefterm(R,x,q,G,tmp,N;l)
    end
    return p
end
function specificFeynmanIntegral(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector, G::graphe ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipV(G, a)
    p=0
    pp=sum(g)
     sz=1
    for k in 1:length(l)
        sz=sz*InvSfunction(z[k],aa)
    end
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * consttermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N)

                elseif f[i][2][j] == 0
                        tmp = tmp * consttermV(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N)

                elseif f[i][2][j] == -2
                        tmp = tmp *looptermV(z[src(ee[j])],q[j],aa,a[j])

                else 
                        tmp = tmp * protermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N)
                end 
        end
            #p=p+f[i][1]*tm
        p=p+f[i][1]*coeftermX(R,x,q,z,G,tmp,N;l)
        if i==length(f)
            pr=sz*p
            p=coefterm2Z(R,x,q,z,pr,g)
        end
    end
    return p
end
function specificFeynmanIntegralo(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector, G::graphe ,a::Vector{Int64} ,o::Vector{Int64};l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipo(G, a,o)
    p=0
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * constterm(x[src(ee[j])], x[dst(ee[j])],N)

                elseif f[i][2][j] == 0
                        tmp = tmp * constterm(x[dst(ee[j])], x[src(ee[j])], N)
                else 
                        tmp = tmp * proterm(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)
                end 
        end
            #p=p+f[i][1]*tm
        p=p+f[i][1]*coefterm(R,x,q,G,tmp,N;l)
    end
    return p
end
function specificFeynmanIntegralo(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector, G::graphe ,a::Vector{Int64} ,o::Vector{Int64};aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipo(G, a,o)
    p=0
    pp=sum(g)
     sz=1
    for k in 1:length(l)
        sz=sz*InvSfunction(z[k],aa)
    end
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * consttermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N)

                elseif f[i][2][j] == 0
                        tmp = tmp * consttermV(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N)

                elseif f[i][2][j] == -2
                        tmp = tmp *looptermV(z[src(ee[j])],q[j],aa,a[j])

                else 
                        tmp = tmp * protermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N)
                end 
        end
            #p=p+f[i][1]*tm
        p=p+f[i][1]*coeftermX(R,x,q,z,G,tmp,N;l)
        
        if i==length(f)
            pr=sz*p
            p=coefterm2Z(R,x,q,z,pr,g)
        end
    end
    return p
end
function feynmanIntegralo(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector, G::graphe,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegralo(R,x,q,G,ai,o;l) 
    end
    return sum
end 
function feynmanIntegralo(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector, G::graphe,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegralo(R,x,q,z,G,ai,o;aa,l,g) 
    end
    return sum
end 
function feynmanIntegral(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector, G::graphe,d::Integer;l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegral(R,x,q,G,ai;l) 
    end
    return sum
end 
function feynmanIntegral(R::Nemo.FmpqMPolyRing, x::Vector, q::Vector,z::Vector, G::graphe,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegral(R,x,q,z,G,ai;aa,l,g) 
    end
    return sum
end 
function feynmanIntegralSum(R::Nemo.FmpqMPolyRing,x::Vector,q::Vector,z::Vector,G::graphe,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    res=0
    for i in 1:d
        res+=feynmanIntegral(R,x,q,z,G,i;aa,l,g)
    end
    return res
end
function feynmanIntegralSum(R::Nemo.FmpqMPolyRing,x::Vector,q::Vector,G::graphe,d::Integer;l=zeros(Int,nv(G)))
    res=0
    for i in 1:d
        res+=feynmanIntegral(R,x,q,G,i;l)
    end
    return res
end
function subt(p::fmpq_mpoly)
coeffs_dict = coeffs(p)
coeffs_array = collect(values(coeffs_dict))
return sum(coeffs_array)
end