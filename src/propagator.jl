@doc raw"""
    constterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, N::Integer)
 
returns the constant term of the propagator

# Examples (without vertex contribution)
 
```julia
julia> constterm(x[1],x[2],3)
 3*x[1]^6 + 2*x[1]^5*x[2] + x[1]^4*x[2]^2
```


    constterm(x1::QQMPolyRingElem, x2::QQMPolyRingElem, z1::QQMPolyRingElem, z2::QQMPolyRingElem,aa::Integer, N::Integer)
    
here `aa=1 ` is the order for sfunction series and `` N=\sum_{n=1}^{3g-3} a_i`` where ``a=[a_1,…,a_{3g-3}]`` is a branch type.

# Examples (without vertex contribution)

```julia
julia> constterm(x[1],x[2],z[1],z[2],1,2)
 1//18*x[1]^4*z[1]^2*z[2]^2 + 1//3*x[1]^4*z[1]^2 + 1//3*x[1]^4*z[2]^2 + 2*x[1]^4 + 1//576*x[1]^3*x[2]*z[1]^2*z[2]^2 + 1//24*x[1]^3*x[2]*z[1]^2 + 1//24*x[1]^3*x[2]*z[2]^2 + x[1]^3*x[2]
```
"""
function constterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, N::Integer)
    p=0
    for i in 1:N
        p=p+i*x1^(N+i)*x2^(N-i)
    end
    return p
end
function proterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer, N::Integer)
    p=0
    for w in 1:a
        if a%w==0
            p = p + w*( x1^( N + w )*x2 ^( N - w ) + x1 ^( N-w )*x2^( N + w ) ) *q^(a)
        end
    end
    return p
end
function constterm(x1::QQMPolyRingElem, x2::QQMPolyRingElem, z1::QQMPolyRingElem, z2::QQMPolyRingElem,aa::Integer, N::Integer)
    p=0
    for i in 1:N
        S1=sfunction(i*z1,aa)
        S2=sfunction(i*z2,aa)
        p=p+S1*S2*i*x1^(N+i)*x2^(N-i)
    end
    return p
end
@doc raw"""
     proterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer, N::Integer)

returns the non constant term of the propagator

# Examples (without vertex contribution)

```julia
julia> proterm(x[1],x[2],q[1],1,2)
x[1]^3*x[2]*q[1] + x[1]*x[2]^3*q[1]
```
# Examples (with vertex contribution)
     proterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem,z1::QQMPolyRingElem, z2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer,aa::Integer, N::Integer)

```julia
julia> proterm(x[1],x[2],z[1],z[2],q[1],1,1,2)
 1//576*x[1]^3*x[2]*q[1]*z[1]^2*z[2]^2 + 1//24*x[1]^3*x[2]*q[1]*z[1]^2 
 + 1//24*x[1]^3*x[2]*q[1]*z[2]^2 + x[1]^3*x[2]*q[1] + 1//576*x[1]*x[2]^3*q[1]*z[1]^2*z[2]^2 
 + 1//24*x[1]*x[2]^3*q[1]*z[1]^2 + 1//24*x[1]*x[2]^3*q[1]*z[2]^2 + x[1]*x[2]^3*q[1]
```
"""
function proterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem,z1::QQMPolyRingElem, z2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer,aa::Integer, N::Integer)
    p=0
    for w in 1:a
        if a%w==0
            S1=sfunction(w*z1,aa)
            S2=sfunction(w*z2,aa)
            p = p + S1*S2*w*( x1^( N + w )*x2 ^( N - w ) + x1 ^( N-w )*x2^( N + w ) ) *q^(a)
        end
    end
    return p
end
