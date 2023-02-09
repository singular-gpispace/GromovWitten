@testset "flip.jl" begin
    ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
    G=graphe(ve)
    R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G))
    a=[1,2,3,2,0,1]
    @test flip(G,a)==[[12, [1, 2, 3, 2, -1, 1]] ,[12, [1, 2, 3, 2, 0, 1]]]
    @test  subt(R,x,q,feynmanIntegral(R,x,q,G ,3))==1792*q[1]^6
    
end