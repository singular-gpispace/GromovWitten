@testset "feynmanint.jl" begin
    ve = [(1, 1), (1, 2), (2, 3), (3, 1)]
    G = FeynmanGraph(ve)
    F = FeynmanIntegral(G)
    q = F.S[3]
    a = [2, 0, 0, 1]
    o = [1, 2, 3]
    g = [1, 0, 0]
    @test feynman_integral_branch_type(F, a, g) == 7 // 4 * q[1]^4 * q[4]^2
    @test feynman_integral_degree_order(F, o, 3, g) == 7 // 8 * q[1]^4 * q[4]^2 + 1 // 8 * q[1]^2 * q[2]^2 * q[3]^2 + 1 // 8 * q[1]^2 * q[2]^2 * q[4]^2 + 1 // 8 * q[1]^2 * q[3]^2 * q[4]^2 + 25 // 8 * q[1]^2 * q[4]^4
    @test feynman_integral_branch_type_order(F, a, o) == 3 * q[1]^4 * q[4]^2
    @test feynman_integral_branch_type(F, a, l=[0, 0, 0]) == 6 * q[1]^4 * q[4]^2
    @test feynman_integral_degree_order(F, o, 3, l=[0, 0, 0]) == 3 * q[1]^4 * q[4]^2 + q[1]^2 * q[2]^2 * q[3]^2 + q[1]^2 * q[2]^2 * q[4]^2 + q[1]^2 * q[3]^2 * q[4]^2 + 9 * q[1]^2 * q[4]^4
    @test sum_of_coeff(feynman_integral_degree_order(F, o, 3)) == 15
    @test sum_of_coeff(feynman_integral_degree(F, 3)) == 72
    v = [(1, 2), (1, 2), (2, 3), (3, 1)]
    G = FeynmanGraph(v)
    F = FeynmanIntegral(G)
    q = F.S[3]
    @test feynman_integral_degree_order(F, o, 3, g) == 5 // 12 * q[1]^4 * q[3]^2 + 11 // 3 * q[1]^4 * q[4]^2 + 5 // 12 * q[1]^2 * q[3]^2 * q[4]^2 + 31 // 4 * q[1]^2 * q[4]^4 + 5 // 12 * q[2]^4 * q[3]^2 + 11 // 3 * q[2]^4 * q[4]^2 + 5 // 12 * q[2]^2 * q[3]^2 * q[4]^2 + 31 // 4 * q[2]^2 * q[4]^4 + 39 // 2 * q[4]^6


end