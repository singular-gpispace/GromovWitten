@testset "feynmanint.jl" begin
    ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
    G=graph(ve)
    R,x,q,z=polynomialring(G,"x","q","z")
    g=[0,0,0]
    a=[2,0,0,1]
    o=[1,2,3]
    @test specificFeynmanIntegralo(x,q,z,G,a,o)==3*q[1]^2*q[4]
    @test specificFeynmanIntegral(x,q,z,G,a,aa=1,l=[0,0,0],g=[0,0,0])==6*q[1]^2*q[4]
    @test feynmanIntegralo(x,q,z,G,o,3)==3*q[1]^2*q[4] + q[1]*q[2]*q[3] + q[1]*q[2]*q[4] + q[1]*q[3]*q[4] + 9*q[1]*q[4]^2
    @test subt(feynmanIntegralo(x,q,z,G,o,3))==15
    @test subt(feynmanIntegral(x,q,z,G ,3))==90
end