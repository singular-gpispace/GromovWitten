@testset "feynmanint.jl" begin
    ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
    G=graph(ve)
    a=[2,0,0,1]
    o=[1,2,3]
    q=Nemo.vars(feynman_integral_degree(G ,3))
    @test sum_of_coeff(feynman_integral_degree(G ,3))==90
    @test feynman_integral_branch_type(G,a,l=[0,0,0])==6*q[1]^4*q[4]^2

    q=Nemo.vars(feynman_integral_degree_order(G,o,3))
    @test feynman_integral_branch_type_order(G,a,o)==3*q[1]^4*q[4]^2
    @test feynman_integral_degree_order(G,o,3,l=[0,0,0])==3*q[1]^4*q[4]^2+ q[1]^2*q[2]^2*q[3]^2+ q[1]^2*q[2]^2*q[4]^2+ q[1]^2*q[3]^2*q[4]^2+ 9*q[1]^2*q[4]^4
    @test sum_of_coeff(feynman_integral_degree_order(G,o,3))==15
  
end