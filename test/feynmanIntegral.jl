@testset "feynmanIntegral.jl" begin
    ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
    G=graphe(ve)
    R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G))
    a=[1,2,3,2,0,1]

    @test specificFeynmanIntegral(R,x,q,G,a)==96*q[1]^2*q[2]^4*q[3]^6*q[4]^4*q[6]^2
    @test feynmanIntegral(R,x,q,G ,2)==8*q[1]^2*q[2]^2 + 8*q[3]^4 + 8*q[4]^4 + 8*q[5]^2*q[6]^2
    @test feynmanIntegralSum(R,x,q,G,2)==8*q[1]^2*q[2]^2 + 8*q[3]^4 + 8*q[4]^4 + 8*q[5]^2*q[6]^2
    @test  subt(R,x,q,feynmanIntegral(R,x,q,G ,3))==1792*q[1]^6
    
end