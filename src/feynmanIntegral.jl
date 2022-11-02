function specificFeynmanIntegral(G::SimpleDiGraph,a::Vector{Int64})
    N = sum(a)
    f = flip(G, a)
   ee=collect(el)
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
               tmp = tmp * term(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)
           end 
       end
      p=p+f[i][1]*coefterm(G,tmp,N)
           
   end
   return p
   
end

function feynmanIntegral(G::SimpleDiGraph,d::Integer)
    a=partition(length(ee),d) 
     sum=0
     for i in 1:length(a)
        sum=sum+specificFeynmanIntegral(G,a[i]) 
     end
     return sum
 end

 function feynmanIntegral(G::SimpleDiGraph,d::Integer)
    a=partition(length(ee),d) 
     sum=0
     for ai in a
        sum=sum+specificFeynmanIntegral(G,ai) 
     end
     return sum
 end

 function feynmanIntegralParallel(G::SimpleDiGraph, d::Integer)
    
    a=partition(d,length(ee)) 
    c=Sys.CPU_THREADS
    while length(a)>0
        res =[]
        for i in 1:c 
            if length(a)>0
                push!(res,specificFeynmanIntegral(G,a[1]))
            a=deleteat!(a,1)
            end
        end
        
    end
    return sum(res)
    
end


function feynmanIntegralSum(G::SimpleDiGraph,d::Integer)
    res=0
   for i in 1:d
        res+=feynmanIntegral(G,i)
    end
    return res
end

function sub(p::fmpq_mpoly)

    if is_constant(R(p))
        return 0
    else
        r=gens(R)
        for i in nv(G)+1:length(r)
           r[i]=q[1]
        end
    end
        return evaluate(p,r)
end