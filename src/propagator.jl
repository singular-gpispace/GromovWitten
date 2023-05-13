

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
function consttermV(x1::fmpq_mpoly, x2::fmpq_mpoly, z1::fmpq_mpoly, z2::fmpq_mpoly,aa::Integer, N::Integer)
    p=0
    for i in 1:N
        S1=Sfunction(i*z1,aa)
        S2=Sfunction(i*z2,aa)
        p=p+S1*S2*i*x1^(N+i)*x2^(N-i)
    end
    return p
end
function protermV( x1::fmpq_mpoly, x2::fmpq_mpoly,z1::fmpq_mpoly, z2::fmpq_mpoly, q::fmpq_mpoly, a::Integer,aa::Integer, N::Integer)
    p=0
    for w in 1:a
        if a%w==0
            S1=Sfunction(w*z1,aa)
            S2=Sfunction(w*z2,aa)
            p = p + S1*S2*w*( x1^( N + w )*x2 ^( N - w ) + x1 ^( N-w )*x2^( N + w ) ) *q^(a)
        end
    end
    return p
end
