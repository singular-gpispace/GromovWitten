@testset "flip.jl" begin
    ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
    G=graph(ve)
    R, x, q = polynomialring(G,"x","q")
    a=[0,2,1,0,0,1]
    @test  subt(feynmanIntegral(x,q,G ,3))==1792
end