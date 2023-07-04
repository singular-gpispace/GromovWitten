function specificFeynmanIntegral( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipV(G, a)
    p=0
    L=lis(G,N,l)
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * constterm(x[src(ee[j])], x[dst(ee[j])],N)

                elseif f[i][2][j] == 0
                        tmp = tmp * constterm(x[dst(ee[j])], x[src(ee[j])], N)
                elseif f[i][2][j] == -2
                        tmp = tmp *looptermV(q[j],a[j])

                else 
                
                        tmp = tmp * proterm(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)

                end 
        end
            #p=p+f[i][1]*tm
       # tmp=filter_term(tmp,x,L)

        p=p+f[i][1]*coeff(tmp,x,L)
       
    end
    return p
end
function specificFeynmanIntegral(x::Vector, q::Vector,z::Vector, G::graph ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipV(G, a)
    L=lis(G,N,l)
    g=2 .* g
    p=0
     sz=1
    for k in 1:length(l)
        sz=sz*filter_term(InvSfunction(z[k],aa),z,g)
    end
    sz=filter_term(sz,z,g)
    for i in 1:length(f)
        tmp=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * filter_term(consttermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N),z,g)

                elseif f[i][2][j] == 0
                        tmp = tmp * filter_term(consttermV(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N),z,g)

                elseif f[i][2][j] == -2
                        tmp = tmp *filter_term(looptermV(z[src(ee[j])],q[j],aa,a[j]),z,g)

                else 
                        tmp = tmp * filter_term(protermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N),z,g)
                end 
        end

    p=p+f[i][1]*coeff(tmp,x,L)
       # p=filter_term(p,z,g)
        if i==length(f)
            pr=sz*p
            p=coeff(pr,z,g)
        end
    end
    return p
end
function specificFeynmanIntegralo( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}  ,o::Vector{Int64};l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipo(G, a,o)
    p=R(0) 
    L=lis(G,N,l)
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * constterm(x[src(ee[j])], x[dst(ee[j])],N)

                elseif f[i][2][j] == 0
                        tmp = tmp * constterm(x[dst(ee[j])], x[src(ee[j])], N)
                elseif f[i][2][j] == -2
                        tmp = tmp *looptermV(q[j],a[j])

                else 
                
                        tmp = tmp * proterm(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)

                end 
        end
            #p=p+f[i][1]*tm
       # tmp=filter_term(tmp,x,L)

        p=p+f[i][1]*coeff(tmp,x,L)
       
    end
    return p
end
function specificFeynmanIntegralo(x::Vector, q::Vector,z::Vector, G::graph ,a::Vector{Int64},o::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = flipo(G, a,o)
    L=lis(G,N,l)
    g=2 .* g
    p=0
     sz=1
    for k in 1:length(l)
        sz=sz*filter_term(InvSfunction(z[k],aa),z,g)
    end
    sz=filter_term(sz,z,g)
    for i in 1:length(f)
        tmp=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * filter_term(consttermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N),z,g)

                elseif f[i][2][j] == 0
                        tmp = tmp * filter_term(consttermV(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N),z,g)

                elseif f[i][2][j] == -2
                        tmp = tmp *filter_term(looptermV(z[src(ee[j])],q[j],aa,a[j]),z,g)

                else 
                        tmp = tmp * filter_term(protermV(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N),z,g)
                end 
        end

    p=p+f[i][1]*coeff(tmp,x,L)
       # p=filter_term(p,z,g)
        if i==length(f)
            pr=sz*p
            p=coeff(pr,z,g)
        end
    end
    return p
end

function feynmanIntegralo( x::Vector, q::Vector, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegralo(x,q,G,ai,o;l) 
    end
    return sum
end 
function feynmanIntegralo( x::Vector, q::Vector,z::Vector, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegralo(x,q,z,G,ai,o;aa,l,g) 
    end
    return sum
end 
function feynmanIntegral( x::Vector, q::Vector, G::graph,d::Integer;l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegral(x,q,G,ai;l) 
    end
    return sum
end 
function feynmanIntegral( x::Vector, q::Vector,z::Vector, G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+specificFeynmanIntegral(x,q,z,G,ai;aa,l,g) 
    end
    return sum
end 
function feynmanIntegralSum(x::Vector,q::Vector,z::Vector,G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    res=0
    for i in 1:d
        res+=feynmanIntegral(x,q,z,G,i;aa,l,g)
    end
    return res
end
function feynmanIntegralSum(x::Vector,q::Vector,G::graph,d::Integer;l=zeros(Int,nv(G)))
    res=0
    for i in 1:d
        res+=feynmanIntegral(x,q,G,i;l)
    end
    return res
end
function subt(p::fmpq_mpoly)
coeffs_dict = coefficients(p)
coeffs_array = collect(values(coeffs_dict))
return sum(coeffs_array)
end