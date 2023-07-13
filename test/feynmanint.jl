@testset "feynmanint.jl" begin
    ve=[(1, 1), (1, 2), (2, 3), (3, 1)]
    G=graph(ve)
    R,x,q,z=polynomial_ring(G,"x","q","z")
    a=[2,0,0,1]
    o=[1,2,3]
    @test feynman_integral_branch_type(x,q,z,G,a,aa=1,g=[1,0,0])==7//4*q[1]^4*q[4]^2
    @test feynman_integral_degree_order(x,q,z,G,o,3,aa=1,g=[1,0,0])==7//8*q[1]^4*q[4]^2 + 1//8*q[1]^2*q[2]^2*q[3]^2 + 1//8*q[1]^2*q[2]^2*q[4]^2 + 1//8*q[1]^2*q[3]^2*q[4]^2 + 25//8*q[1]^2*q[4]^4
    @test feynman_integral_branch_type_order(x,q,G,a,o)==3*q[1]^4*q[4]^2
    @test feynman_integral_branch_type(x,q,G,a,l=[0,0,0])==6*q[1]^4*q[4]^2
    @test feynman_integral_degree_order(x,q,G,o,3,l=[0,0,0])==3*q[1]^4*q[4]^2+ q[1]^2*q[2]^2*q[3]^2+ q[1]^2*q[2]^2*q[4]^2+ q[1]^2*q[3]^2*q[4]^2+ 9*q[1]^2*q[4]^4
    @test sum_of_coeff(feynman_integral_degree_order(x,q,G,o,3))==15
    @test sum_of_coeff(feynman_integral_degree(x,q,G ,3))==90
    ve=[]
    G=graph(ve)
    @test_throws DomainError R,x,q=polynomial_ring(G,"x","q")


end