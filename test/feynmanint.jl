@testset "feynmanint.jl" begin
    ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
    G=graph(ve)
    R,x,q=polynomial_ring(G,"x","q")
    a=[2,0,0,1]
    o=[1,2,3]
    @test loopterm(q[1],2)==3*q[1]^4
    @test feynman_integral_branch_type_order(x,q,G,a,o)==3*q[1]^4*q[4]^2
    @test feynman_integral_branch_type(x,q,G,a,l=[0,0,0])==6*q[1]^4*q[4]^2
    @test feynman_integral_degree_order(x,q,G,o,3,l=[0,0,0])==3*q[1]^4*q[4]^2+ q[1]^2*q[2]^2*q[3]^2+ q[1]^2*q[2]^2*q[4]^2+ q[1]^2*q[3]^2*q[4]^2+ 9*q[1]^2*q[4]^4
    @test sum_of_coeff(feynman_integral_degree_order(x,q,G,o,3))==15
    @test sum_of_coeff(feynman_integral_degree(x,q,G ,3))==90

end