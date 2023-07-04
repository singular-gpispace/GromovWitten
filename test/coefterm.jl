@testset "flip.jl" begin
    ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
    G=graph(ve)
    R, x, q = polynomial_ring(G,"x","q")
    a=[0,2,1,0,0,1]
    @test  subt(feynman_integral_degree(x,q,G ,3))==1792
end