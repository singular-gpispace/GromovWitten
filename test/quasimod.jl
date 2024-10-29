@testset "quasimod.jl" begin
    w = [2, 4, 6]
    d = 12
    @test number_of_monomials( d) == 22
    import Nemo: QQFieldElem
    S, (E2, E4, E6) = polynomial_ring(QQ, ["E2", "E4", "E6"])
    R, q = polynomial_ring(QQ, ["q"])
    ve = [(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)]
    G = FeynmanGraph(ve)
    F = FeynmanIntegral(G)
    @testset "express_as_eisenstein_series test" begin
        expected_result =Any[E6, E6^2, E4, E4*E6, E4^2, E4^3, E2, E2*E6, E2*E4, E2*E4*E6, E2*E4^2, E2^2, E2^2*E6, E2^2*E4, E2^2*E4^2, E2^3, E2^3*E6, E2^3*E4, E2^4, E2^4*E4, E2^5, E2^6]
        @test express_as_eisenstein_series(12) == expected_result
    end

    @test sum_of_divisor_powers(1, 1) == 1
    @testset "sum_of_divisor_powers test" begin
        @test_throws ErrorException sum_of_divisor_powers(0, 1)
    end
    @testset "eisenstein_series even" begin
        @test_throws ErrorException eisenstein_series(6, 3)
    end

    @test eisenstein_series(6, 2) == -288 * q[1]^12 - 144 * q[1]^10 - 168 * q[1]^8 - 96 * q[1]^6 - 72 * q[1]^4 - 24 * q[1]^2 + 1


    @testset "express_as_powers test" begin
        expected_result = [
            60480 * q[1]^12 + 30240 * q[1]^10 + 17520 * q[1]^8 + 6720 * q[1]^6 + 2160 * q[1]^4 + 240 * q[1]^2 + 1
            -288*q[1]^12 - 144*q[1]^10 - 168*q[1]^8 - 96*q[1]^6 - 72*q[1]^4 - 24*q[1]^2 + 1
            82944 * q[1]^24 + 82944 * q[1]^22 + 117504 * q[1]^20 + 103680 * q[1]^18 + 97344 * q[1]^16 + 66816 * q[1]^14 + 39744 * q[1]^12 + 21600 * q[1]^10 + 9456 * q[1]^8 + 3264 * q[1]^6 + 432 * q[1]^4 - 48 * q[1]^2 + 1
        ]
        @test express_as_powers(6, 4) == expected_result
    end

    ep = express_as_powers(6, 4)
    @testset "filter_vector test" begin
        expected_result = [
            6720 * q[1]^6 + 2160 * q[1]^4 + 240 * q[1]^2 + 1
            -96*q[1]^6 - 72*q[1]^4 - 24*q[1]^2 + 1
            3264 * q[1]^6 + 432 * q[1]^4 - 48 * q[1]^2 + 1
        ]
        @test filter_vector(ep, q, 6) == expected_result
    end

    @testset "polynomial_to_matrix test" begin
        eq = filter_vector(ep, q, [2])
        expected_result = [
            1 1 1
            0 0 0
            240 -24 -48]
        result = polynomial_to_matrix(eq)

        for idx in Tuple(eachindex(result))
            i, j = Tuple(idx)
            @test result[i, j] == QQFieldElem(expected_result[i, j])
        end
    end

    @testset "matrix_of_integral " begin
        Iq3 = substitute(feynman_integral_degree_sum(F, 3))
        expected_result = [
            0
            0
            0
            0
            32
            0
            1792
        ]
        result = matrix_of_integral(Iq3)

        for idx in Tuple(eachindex(result))
            i, j = Tuple(idx)
            @test result[i, j] == QQFieldElem(expected_result[i, j])
        end
    end

    @testset "solve_polynomial_system Error" begin
        max_degree = 8
        d = 4
        weightmax = 12
        A = polynomial_to_matrix(filter_vector(express_as_powers(max_degree, 4), q, [max_degree]))
        Q = matrix_of_integral(substitute(feynman_integral_degree_sum(F, d)))
        expected_result = "The system has no solution"

        result = solve_polynomial_system(A, Q, side=:right)
        @test result == expected_result
    end

    @testset "solve_polynomial_system" begin
        max_degree = 46
        d = 6
        weightmax = 12
        A = polynomial_to_matrix(filter_vector(express_as_powers(max_degree, weightmax), q, [max_degree]))
        Iq=61425005056*q[1]^46+43646419584*q[1]^44+29331341312*q[1]^42+20067375616*q[1]^40 + 12961886976*q[1]^38 + 8490271392*q[1]^36 + 5225373696*q[1]^34 + 3233267712*q[1]^32 + 1875116544*q[1]^30 + 1079026432*q[1]^28 + 577972224*q[1]^26 + 302347264*q[1]^24 + 145337600*q[1]^22 + 66497472*q[1]^20 + 27353088*q[1]^18 + 10246144*q[1]^16 + 3294720*q[1]^14 + 886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
        Q = matrix_of_integral(Iq)
        expected_result = (
            QQFieldElem(1 // 93312),
            QQFieldElem[0; 4; 0; 0; 0; 4; 0; 0; 0; -12; 0; 0; 0; 0; -3; 0; 4; 0; 0; 6; 0; -3;;]
        )

        result = solve_polynomial_system(A, Q)
        @test result[1] == expected_result[1]

        for i in size(expected_result[2], 1)
            @test result[2][i] == expected_result[2][i]
        end
    end
 
    @testset "quasimodular_form" begin
        expected_result = (1//93312, -3*E2^6 + 6*E2^4*E4 + 4*E2^3*E6 - 3*E2^2*E4^2 - 12*E2*E4*E6 + 4*E4^3 + 4*E6^2)
        ve = [(1, 2), (1, 2), (1, 2)]

        F = FeynmanIntegral(ve)
        weightmax=12
        m=number_of_monomials(6)
        Iq=61425005056*q[1]^46+43646419584*q[1]^44+29331341312*q[1]^42+20067375616*q[1]^40 + 12961886976*q[1]^38 + 8490271392*q[1]^36 + 5225373696*q[1]^34 + 3233267712*q[1]^32 + 1875116544*q[1]^30 + 1079026432*q[1]^28 + 577972224*q[1]^26 + 302347264*q[1]^24 + 145337600*q[1]^22 + 66497472*q[1]^20 + 27353088*q[1]^18 + 10246144*q[1]^16 + 3294720*q[1]^14 + 886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
        result=quasimodular_form(Iq,weightmax)
        @test result == expected_result

    end
    @testset "quasimodularity_form" begin
        expected_result = (1//2160, 5*E2^3 - 3*E2*E4 - 2*E6)
        ve = [(1, 2), (1, 2), (1, 2)]

        F = FeynmanIntegral(ve)
        weightmax=6
        m=number_of_monomials(6)
        Iq=substitute(feynman_integral_degree_sum(F, m))
        result=quasimodularity_form(Iq,weightmax)
        @test result == expected_result

    end

end