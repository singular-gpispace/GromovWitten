@testset "feynmanIntegral.jl" begin
    ve=[(1, 3), (1,2),(1, 2), (2, 4) ,(3, 4),(3,4)]
    G=graph(ve)
    R,x,q=polynomialring(G,"x","q")
    a=[0,2,1,0,0,1]
    o=[1,3,4,2]
    @test specificFeynmanIntegralo(R,x,q,G,a,o)==128*q[2]^2*q[3]*q[6]
    @test specificFeynmanIntegral(R,x,q,G,a)==256*q[2]^2*q[3]*q[6]
    @test feynmanIntegral(R,x,q,G ,2)==8*q[1]^2 + 8*q[2]*q[3] + 8*q[4]^2 + 8*q[5]*q[6]
    @test subt( feynmanIntegralo(R,x,q,G,o,4))==3088
    @test subt(feynmanIntegral(R,x,q,G,4))==25344
end