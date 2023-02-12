
function coefterm(G::graphe ,p::fmpq_mpoly,d::Integer)

    ee=Edge.(G.edge)
    G=DiGraph(Edge.(G.edge))

    L=zeros(Int,nv(G))
    for ev in ee
        L[src(ev)]=L[src(ev)]+d
        L[dst(ev)]=L[dst(ev)]+d
    end
    for (i, l) in enumerate(L)         
        f=p
        m=coefficients(f,i)
        if size(m,1)<=l
            return 0 
        end
        p=m[l+1,1]
    end
    
    return p
end 

function partition(n::Integer, k::Integer)
    if k == 1
        return [[n]]
    end
    if n == 0
        return [[0] * k]
    end
    result = Array{Int64}[]
    for A in with_replacement_combinations(1:n, k)
        item = zeros(Int64, n)
        for a in A
            item[a] += 1
        end
        push!(result, item)
    end
    return result
end
function preimg(S::Vector, a::Integer)
    for (i,s) in enumerate(S)
        if s==a
            return i
            break
        end
    end
end

function sgn(G::graphe ,p::Vector,a::Vector) #graph G, list p, branch type a
    
    ee=Edge.(G.edge)

    b=zeros(Int,length(a))
    for (i,(ai,ev)) in enumerate(zip(a,ee))
       if ai==0
           if preimg(p,src(ev))<preimg(p,dst(ev))
                    b[i]=-1
                else 
                    b[i]=0
            end
            else
                b[i]=ai
        end         
    end
            return b
end


function flip( G::graphe,  a::Vector)

    ee=Edge.(G.edge)

    p=[]
    b=Vector{Vector{Any}}()

    l=zeros(Int,nv(G))
    y=Vector{Any}()

    for (ev,ai) in zip(ee,a)
        if ai==0
            l[src(ev)] =1
            l[dst(ev)] =1
        end
    end

    for (i,li) in enumerate(l)
           if li==1
            push!(p,i)
        end

    end   
    p=collect(permutations(p))

    for ga in p 
        push!(y,sgn(G,ga,a))    
    end
    dd=div(factorial( nv(G) ) , length( p ) )
 
    for (key, val) in countmap(y)
        push!(b,[ val*dd, key]) 
    end
    return b   
end