#=
"""
  derivative_series_from_initial!(; initial_series::ComplexOffsetVector, 
    derivatives_series ::ComplexOffsetMatrix)

Transform initial series into array of derivatives' and integrals' series

Given long Taylor series
``f(z) = \\sum_{n=0}^N a_n z^n``
generates array of series
``\\frac{1}{k!}f^{(k)}(z) = \\sum_{n=0}^M b_{k,n} z^n, 0 < z ≤ K_2``,
``f^{(k)}(z) = \\sum_{n=0}^M b_{k,n} z^n, K_1 ≤ z < 0``,
where K1 and K2 are first and last indices if *derivative_series*.

Assuming  N ≥ M + (K_2 - K_1 + 1).
"""
function derivative_series_from_initial!(;
        initial_series     ::ComplexOffsetVector,
        derivatives_series ::ComplexOffsetMatrix
        )
    #####
    upper_term       = lastindex(derivatives_series, 1)
    upper_derivative = lastindex(derivatives_series, 2)
    lower_derivative = firstindex(derivatives_series, 2)

    #derivatives_series[0 : upper_term, 0] .= initial_series[0 : upper_term]
    for n in 0 : upper_term
        derivatives_series[n, 0] = initial_series[n]
    end

    for k in 1 : upper_derivative
        for n in 1 : (upper_term + upper_derivative - k + 1)
            initial_series[n-1] = initial_series[n] * n / (k+1)
        end
        initial_series[upper_term + upper_derivative - k + 1] = .0im

        #@views derivatives_series[:, k] .= initial_series[0 : upper_term]
        for n in 0 : upper_term
            derivatives_series[n, k] = initial_series[n]
        end
    end

    for k in -1 : -1 : lower_derivative
        derivatives_series[0, k] = initial_series[k]
        for n in 1 : upper_term
            derivatives_series[n, k] = derivatives_series[n-1, k+1] / n
        end
    end
end

function construct_dict_of_wei_series(;
    delta_to_nterms ::RComplex2IntDict,
    wei             ::Weierstrass
    )
#####
    max_n_terms = 0
    for v in values(delta_to_nterms)
        if v[1] > max_n_terms
            max_n_terms = v[1]
        end
    end

    wei_sequence = OffsetArray{ComplexF64}(undef, -2 : 2*max_n_terms)

    wei_series_dict = Dict{RationalComplex, ComplexOffsetMatrix}()

    for delta in keys(delta_to_nterms)
        upper_term = delta_to_nterms[delta]

        # Для вычисления разложений в ряд производных нужно разложить исходную
        # функцию с запасом членов.
        # Всего нужно n производных по n членов в каждой.
        diff_upper_term = 2 * upper_term

        if delta != 0
            z = real(delta)*wei.ω1 + imag(delta)*wei.ω3
            norm_factor = abs(z)

            weierstrass_normalized!(z=z,
                                norm_factor=norm_factor,
                                w=wei,
                                upper_term=diff_upper_term,
                                output=wei_sequence)
        else
            expand_at_pole!(output=wei_sequence,
                            weierstrass_data=wei,
                            upper_term=diff_upper_term)
        end

        series_array = OffsetArray{ComplexF64}(undef, 0:upper_term, -1:upper_term-2)

        derivative_series_from_initial!(initial_series=wei_sequence,
                                        derivatives_series=series_array)

        wei_series_dict[delta] = series_array
    end

    return wei_series_dict
end
=#