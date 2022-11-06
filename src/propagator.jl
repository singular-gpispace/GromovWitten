

function constterm( x1::fmpq_mpoly, x2::fmpq_mpoly, N::Integer)
    p=0
    for i in 1:N
        p=p+i*x1^(N+i)*x2^(N-i)
    end
    return p
end

function proterm( x1::fmpq_mpoly, x2::fmpq_mpoly, q::fmpq_mpoly, a::Integer, N::Integer)
    p=0
    for w in 1:a
        if a%w==0
            p = p + w*( x1^( N + w )*x2 ^( N - w ) + x1 ^( N-w )*x2^( N + w ) ) *q^(2*a)

        end
    end
    return p
end

function propagator( x1::fmpq_mpoly, x2::fmpq_mpoly,q::fmpq_mpoly, N::Integer)
    p=constterm(x1,x2,N)
    for i in 1:N
        p=p+proterm(x1,x2,q,i,N)
    end
    return p
end

