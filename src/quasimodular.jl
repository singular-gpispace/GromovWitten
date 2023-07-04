function express_as_eisenstein_series(max_degree::Int)
    expressions = []
    for e2 in 0:max_degree
        for e4 in 0:max_degree
            for e6 in 0:max_degree
                degree = 2*e2 + 4*e4 + 6*e6
                if degree == max_degree
                    push!(expressions, :(E2^$e2 * E4^$e4 * E6^$e6))
                end
            end
        end
    end
    result = Expr(:call, :+, expressions...)
    return expressions
end
function filter_term(pols::Union{QQMPolyRingElem, Int64}, variables::Vector{QQMPolyRingElem}, power::Vector{Int64})
    T = parent(variables[1])
    if typeof(pols) <: Integer
        pols = T(pols)
    end
    gensR = gens(T)
    position = [findfirst(var -> var == vi, gensR) for vi in variables]
    result = zero(pols)
    d = Vector{Int}(undef, length(variables))  # Preallocate d
    @inbounds for term in terms(T(pols))
        for j in 1:length(variables)
            po = position[j]
            d[j] = degree_fmpz(term, po)
        end
        if all(d .<= power)
            result += term
        end
    end
    
    return result
end
function filter_vector(polyvector::Vector{QQMPolyRingElem}, variables::Vector{QQMPolyRingElem}, power::Vector{Int64})
    result = Vector{QQMPolyRingElem}()
    for pols in polyvector
        push!(result,filter_term(pols, variables::Vector{QQMPolyRingElem}, power::Vector{Int64}))
    end
    return result
end
function sum_of_divisor_powers(n::Int, k::Int)
    divs = divisors(n)
    div_sum = sum(div^k for div in divs)
    return div_sum
end
# Compute the Eisenstein series E2(q[1]) up to a specified number of terms
function eisenstein_series(q, num_terms,k)   
    e2 = 1 - ((((2*k)// bernoulli(k))) )* sum(sum_of_divisor_powers(d, k-1) * q^(2*d) for d in 1:num_terms)
    return e2
end  

function express_as_powers(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}}, max_degree)
    if max_degree % 2 != 0
        error("input must be even")
    end
    result = Vector{QQMPolyRingElem}()
    nb=Int64(max_degree//2)
    E6=eisenstein_series(q, nb,6)
    E4=eisenstein_series(q, nb,4)
    E2=eisenstein_series(q, nb,2)
    for e2 in 0:max_degree
        for e4 in 0:max_degree
            for e6 in 0:max_degree
                degree = 2*e2 + 4*e4 + 6*e6
                if degree == max_degree
                    push!(result,  (E2^e2) * (E4^e4) * (E6^e6))
                end
            end
        end
    end
    return result
end
function polynomial_to_matrix(vect::Vector{QQMPolyRingElem})
    A = matrix(QQ, hcat([collect(coefficients_of_univariate(poly)) for poly in vect]...))
    return A
end
function matrix_of_integral(Iq::fmpq_mpoly) 
# Obtain the coefficients of Iq
Q = coefficients_of_univariate(Iq)
    
# Convert Q to a matrix representation
Q_matrix = Matrix(transpose(reduce(hcat, Q)))

# Convert the matrix to the QQ type so we can call can_with_solution.
Q = matrix(QQ, Q_matrix)
    return Q
end
function solve_polynomial_system(A::QQMatrix, Q::QQMatrix)
    m, n = size(A)
    
    # Solve the polynomial system
    can_solve, x = can_solve_with_solution(A, Q)
    
    # Compute the common factor and scaled vector
    if can_solve
        common_factor = gcd(x...)
        coeff_vector = [val // common_factor for val in x]
        return  common_factor, coeff_vector
    else
        return can_solve,x
    end
end
    function quasi_matrix(q::Union{QQMPolyRingElem, Vector{QQMPolyRingElem}},Iq::QQMPolyRingElem, max_degree::Int64)
    if typeof(q) == Vector{QQMPolyRingElem}
        Evector=filter_vector(express_as_powers( q[1],max_degree),[q[1]],[max_degree])
        A=polynomial_to_matrix(Evector)
        Q=matrix_of_integral(Iq)
            return solve_polynomial_system(A,Q)
    else
        Evector=filter_vector(express_as_powers( q,max_degree),[q],[max_degree])
        A=polynomial_to_matrix(Evector)
        Q=matrix_of_integral(Iq)
            return solve_polynomial_system(A,Q)
    end
end