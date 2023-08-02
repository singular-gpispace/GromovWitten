@testset "quasimod.jl" begin
    import Nemo: QQFieldElem

    ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
    G=graph(ve)
    q=vars(feynman_integral_degree_sum(G,6))
    @test substitute(feynman_integral_degree_sum(G,6))==886656*q[1]^12 + 182272*q[1]^10 + 25344*q[1]^8 + 1792*q[1]^6 + 32*q[1]^4
    Iq=substitute(feynman_integral_degree_sum(G,6))
   # R, x, q = polynomial_ring(G,"x","q")

    S,E2,E4,E6=@polynomial_ring(QQ,E2,E4,E6)


    @testset "express_as_eisenstein_series test" begin
        expected_result = [
            E6^2
            E4^3
            E2*E4*E6
            E2^2*E4^2
            E2^3*E6
            E2^4*E4
            E2^6
        ]
        
        @test express_as_eisenstein_series(E2,E4,E6,12) == expected_result
    end
    
    @test sum_of_divisor_powers(1,1)==1
    @testset "sum_of_divisor_powers test" begin
        @test_throws ErrorException sum_of_divisor_powers(0, 1)
    end    
    @testset "eisenstein_series even" begin
        @test_throws ErrorException eisenstein_series(q[1],6,3)
    end  

    @test eisenstein_series(q[1],6,2)==-288*q[1]^12 - 144*q[1]^10 - 168*q[1]^8 - 96*q[1]^6 - 72*q[1]^4 - 24*q[1]^2 + 1
    
    @testset "express_as_powers even" begin
        @test_throws ErrorException express_as_powers(q,3)
    end  

    @testset "express_as_powers test" begin
        expected_result = [
            -122976*q[1]^6 - 16632*q[1]^4 - 504*q[1]^2 + 1
            -645120*q[1]^12 - 691200*q[1]^10 - 339840*q[1]^8 - 62496*q[1]^6 - 3672*q[1]^4 + 216*q[1]^2 + 1
            -884736*q[1]^18 - 1990656*q[1]^16 - 2156544*q[1]^14 - 1340928*q[1]^12 - 497664*q[1]^10 - 95040*q[1]^8 - 3744*q[1]^6 + 1512*q[1]^4 - 72*q[1]^2 + 1
            
        ]
        
        @test  express_as_powers(q,6) == expected_result
    end

     ep=express_as_powers(q,6)
     @testset "filter_vector test" begin
        expected_result = [
            -122976*q[1]^6 - 16632*q[1]^4 - 504*q[1]^2 + 1
            -62496*q[1]^6 - 3672*q[1]^4 + 216*q[1]^2 + 1
            -3744*q[1]^6 + 1512*q[1]^4 - 72*q[1]^2 + 1
        ]
        @test filter_vector(ep,q,6)==expected_result
    end
    
    @testset "polynomial_to_matrix test" begin
        eq = filter_vector(ep, q, [2])
        expected_result = [
            1     1     1
            0     0     0
         -504   216   -72
         
        ]
        result = polynomial_to_matrix(eq)
        
        for idx in Tuple(eachindex(result))
            i, j = Tuple(idx)
            @test result[i, j] == QQFieldElem(expected_result[i, j])
        end
    end
    
    @testset "matrix_of_integral " begin
        Iq3=substitute(feynman_integral_degree_sum(G,3))
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
        max_degree=8
        d=4
        A = polynomial_to_matrix(filter_vector(express_as_powers(q, max_degree), q, [max_degree]))
        Q = matrix_of_integral(substitute( feynman_integral_degree_sum( G, d)))
        expected_result = "The system has no solution"
    
        result = solve_polynomial_system(A, Q)
        @test result == expected_result
    end

    @testset "solve_polynomial_system" begin
        max_degree=12
        d=6
        A = polynomial_to_matrix(filter_vector(express_as_powers(q, max_degree), q, [max_degree]))
        Q = matrix_of_integral(substitute( feynman_integral_degree_sum( G, d)))
        expected_result = (QQFieldElem(1//93312),QQFieldElem[4
        4
        -12
        -3
        4
        6
        -3] )
    
        result = solve_polynomial_system(A, Q)
        @test result[1] == expected_result[1]

        for i in size(expected_result[2],1)
            @test result[2][i] == expected_result[2][i]
        end
    end

    @testset "quasi_matrix case QQFieldElem " begin
        expected_result = (QQFieldElem(1//93312),QQFieldElem[4
        4
        -12
        -3
        4
        6
        -3] )
    
        result = quasi_matrix(q[1],Iq,12)
        @test result[1] == expected_result[1]

        for i in size(expected_result[2],1)
            @test result[2][i] == expected_result[2][i]
        end
    end
    
    @testset "quasi_matrix case Vector{QQMPolyRingElem}" begin
        expected_result = (QQFieldElem(1//93312),QQFieldElem[4
        4
        -12
        -3
        4
        6
        -3] )
    
        result = quasi_matrix(q,Iq,12)
        @test result[1] == expected_result[1]

        for i in size(expected_result[2],1)
            @test result[2][i] == expected_result[2][i]
        end
    end
    Jq= 843264*q[1]^12 + 165888*q[1]^10 + 20736*q[1]^8 + 1152*q[1]^6
    @testset "quasi_matrix error" begin
        @test_throws DimensionMismatch quasi_matrix(q,Jq,8)
    end  
    @testset "quasimodular_form" begin
        expected_result = (1//20736, -E2^6 + 3*E2^4*E4 - 3*E2^2*E4^2 + E4^3)
    
        result = quasimodular_form(E2,E4,E6,q,Jq,12)
        @test result == expected_result
    
    end

end