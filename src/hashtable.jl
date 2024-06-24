function cache_integral_result(F::FeynmanIntegral, integral_type::Symbol, args::Vector{Int64}, result::QQMPolyRingElem)
    if !haskey(F.integral_cache, integral_type)
        F.integral_cache[integral_type] = Dict{Vector{Int64},QQMPolyRingElem}()
    end
    F.integral_cache[integral_type][args] = result
end
function get_integral_from_cache(F::FeynmanIntegral, integral_type::Symbol, args::Vector{Int64})
    if haskey(F.integral_cache, integral_type)
        integral_dict = F.integral_cache[integral_type]
        if haskey(integral_dict, args)
            return integral_dict[args]
        end
    end
    return nothing
end

# Function to compute and cache the branch type integral
function feynman_integral_branch_type_cache(F::FeynmanIntegral, a::Vector{Int64}; l=zeros(Int, nv(F.G)))
    # Check if the integral has been computed before
    cached_result = get_integral_from_cache(F, :branch, a)
    if cached_result !== nothing
        return cached_result
    end

    # If not computed before, compute the integral
    result = feynman_integral_branch_type(F, a; l)

    # Cache the result for future use
    cache_integral_result(F, :branch, a, result)

    return result
end

# Function to compute and cache the degree integral
function feynman_integral_degree_cache(F::FeynmanIntegral, d::Int64; l=zeros(Int, nv(F.G)))
    # Check if the integral has been computed before
    cached_result = get_integral_from_cache(F, :degree, [d])
    if cached_result !== nothing
        return cached_result
    end

    # If not computed before, compute the integral
    result = feynman_integral_degree(F, d; l)

    # Cache the result for future use
    cache_integral_result(F, :degree, [d], result)

    return result
end
function feynman_integral_degree_sum_cache(F::FeynmanIntegral, d::Union{Int64,Vector{Int64}}; l=zeros(Int, nv(F.G)))
    # Check if the current degree is already computed and get the cached result
    current_result = get_integral_from_cache(F, :sum, [d])
    if current_result !== nothing
        return current_result
    else
        # If the current degree is not cached, find the maximal b < d with cached result
        b = 0
        for i in reverse(1:d[1]-1)
            cached_result = get_integral_from_cache(F, :sum, [i])
            if cached_result !== nothing
                b = i
                break
            end
        end

        # Compute the sum for the current degree [b+1, d]
        sum_result = substitute(feynman_integral_degree_sum(F, [b + 1, d]))
        # Get the cached result for the maximal b
        cached_b_result = get_integral_from_cache(F, :sum, [b])

        # Calculate the result for the current degree [b, d]
        if cached_b_result !== nothing
            current_result = sum_result + cached_b_result
        else
            current_result = sum_result
        end

        # Cache the result for the current degree
        cache_integral_result(F, :sum, [d], current_result)
    end

    return current_result
end