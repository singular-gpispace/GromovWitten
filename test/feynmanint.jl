@testset "feynmanint.jl" begin
    ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
    G=graphe(ve)
    R,x,q,z=polynomialring(G,"x","q","z")
    g=[0,0,0]
    a=[2,0,0,1]
    o=[1,2,3]
    @test specificFeynmanIntegralo(R,x,q,z,G,a,o)==3*q[1]^2*q[4]
    @test specificFeynmanIntegral(R,x,q,z,G,a,aa=1,l=[0,0,0],g=[0,0,0])==6*q[1]^2*q[4]
    @test feynmanIntegralo(R,x,q,z,G,o,3)==3*q[1]^2*q[4] + q[1]*q[2]*q[3] + q[1]*q[2]*q[4] + q[1]*q[3]*q[4] + 9*q[1]*q[4]^2
    @test subt(feynmanIntegralo(R,x,q,z,G,o,3))==15
    @test subt(feynmanIntegral(R,x,q,z,G ,3))==90
end