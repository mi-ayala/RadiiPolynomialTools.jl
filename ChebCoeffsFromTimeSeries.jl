using ApproxFun


"""
    get_cheb_coefficients(m)

Calculate the Chebyshev coefficients for a given `m`.

"""
function get_cheb_coefficients(t, f_t, m)
    S = ApproxFun.Chebyshev(ApproxFun.Interval(-1, 1))
    n = length(t)

    V = Array{Float64}(undef, n, m)
    for k = 1:m
        V[:, k] = Fun(S, [zeros(k - 1); 1]).(t)
    end
    Fun(S, V \ f_t)

end

"""
    find_best_m(start_m, end_m, step_m, threshold)

Find the `m` value that results in the smallest decay when passed to `get_cheb_coefficients`.

"""
function find_best_Cheb_mode(t, f_t, m_max)

    m = 5
    decay = ApproxFun.coefficients(get_cheb_coefficients(t, f_t, m))[end]
    println("m", " ", "Decay")
    while abs(decay) > 1e-12 && m < m_max
        m += 1
        decay = ApproxFun.coefficients(get_cheb_coefficients(t, f_t, m))[end]
        println(m, " ", decay)
    end

    return m

end

function cheb_2_scaling(s)

    s = [s[1]; (1 / 2) * s[2:end]]

end

function cheb_coeffs_from_time_series(t, f_t, m_max)

    m = find_best_Cheb_mode(t, f_t, m_max)
    cheb_2_scaling(ApproxFun.coefficients(get_cheb_coefficients(t, f_t, m)))

end

function get_cheb_sequence(sol, max_M)

    s1 = cheb_coeffs_from_time_series(sol.t, sol[1, :], max_M)
    s2 = cheb_coeffs_from_time_series(sol.t, sol[2, :], max_M)
    s3 = cheb_coeffs_from_time_series(sol.t, sol[3, :], max_M)
    s4 = cheb_coeffs_from_time_series(sol.t, sol[4, :], max_M)
    s5 = cheb_coeffs_from_time_series(sol.t, sol[5, :], max_M)
    s6 = cheb_coeffs_from_time_series(sol.t, sol[6, :], max_M)
    s7 = cheb_coeffs_from_time_series(sol.t, sol[7, :], max_M)
    s8 = cheb_coeffs_from_time_series(sol.t, sol[8, :], max_M)

    sizes = [length(s1) length(s2) length(s3) length(s4) length(s5) length(s6) length(s7) length(s8)]
    N_C = (maximum(sizes) - 1) + 3

    solitonSpace = RadiiPolynomial.Chebyshev(N_C)^8
    s = Sequence(solitonSpace, zeros(Float64, dimension(solitonSpace)))

    RadiiPolynomial.component(s, 1)[0:sizes[1]-1] .= s1
    RadiiPolynomial.component(s, 2)[0:sizes[2]-1] .= s2
    RadiiPolynomial.component(s, 3)[0:sizes[3]-1] .= s3
    RadiiPolynomial.component(s, 4)[0:sizes[4]-1] .= s4
    RadiiPolynomial.component(s, 5)[0:sizes[5]-1] .= s5
    RadiiPolynomial.component(s, 6)[0:sizes[6]-1] .= s6
    RadiiPolynomial.component(s, 7)[0:sizes[7]-1] .= s7
    RadiiPolynomial.component(s, 8)[0:sizes[8]-1] .= s8

    return s


end



##### get_cheb_sequence(sol, max_M) where sol is the output of the integration with DifferentialEquations.jl
