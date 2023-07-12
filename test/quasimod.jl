@testset "quasimod.jl" begin
    ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
    G=graph(ve)
    R, x, q = polynomial_ring(G,"x","q")
    @test substitute(q,feynman_integral_degree_sum(x,q,G,6))==886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
    f=substitute(q,feynman_integral_degree_sum(x,q,G,6))

    @testset "express_as_eisenstein_series test" begin
        expected_result = [
            :(E2 ^ 0 * E4 ^ 0 * E6 ^ 2),
            :(E2 ^ 0 * E4 ^ 3 * E6 ^ 0),
            :(E2 ^ 1 * E4 ^ 1 * E6 ^ 1),
            :(E2 ^ 2 * E4 ^ 2 * E6 ^ 0),
            :(E2 ^ 3 * E4 ^ 0 * E6 ^ 1),
            :(E2 ^ 4 * E4 ^ 1 * E6 ^ 0),
            :(E2 ^ 6 * E4 ^ 0 * E6 ^ 0)
        ]
        
        @test express_as_eisenstein_series(12) == expected_result
    end
    
    @test sum_of_divisor_powers(1,1)==1
    @test eisenstein_series(q[1],6,2)==-288*q[1]^12 - 144*q[1]^10 - 168*q[1]^8 - 96*q[1]^6 - 72*q[1]^4 - 24*q[1]^2 + 1
    @testset "express_as_powers test" begin
        expected_result = [
            -122976*q[1]^6 - 16632*q[1]^4 - 504*q[1]^2 + 1
            -645120*q[1]^12 - 691200*q[1]^10 - 339840*q[1]^8 - 62496*q[1]^6 - 3672*q[1]^4 + 216*q[1]^2 + 1
            -884736*q[1]^18 - 1990656*q[1]^16 - 2156544*q[1]^14 - 1340928*q[1]^12 - 497664*q[1]^10 - 95040*q[1]^8 - 3744*q[1]^6 + 1512*q[1]^4 - 72*q[1]^2 + 1
            
        ]
        
        @test  express_as_powers(q,6) == expected_result
    end

    ep=express_as_powers(q,6)
   #= @testset "filter_vector test" begin
        expected_result = [
            -122976*q[1]^6 - 16632*q[1]^4 - 504*q[1]^2 + 1
            -62496*q[1]^6 - 3672*q[1]^4 + 216*q[1]^2 + 1
            -3744*q[1]^6 + 1512*q[1]^4 - 72*q[1]^2 + 1
        ]
        @test filter_vector(ep,q,[6])==expected_result
    end
    @testset "polynomial_to_matrix test" begin
        expected_result = [
            [1 1 1; 0 0 0; -504 216 -72; 0 0 0; -16632 -3672 1512; 0 0 0; -122976 -62496 -3744]
        ]
        @test polynomial_to_matrix( filter_vector(ep,q,[6]))==expected_result
    end =#
end