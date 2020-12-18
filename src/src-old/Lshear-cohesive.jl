#=
using OffsetArrays

using ..Ellipticals: Weierstrass, term_expansion_normalized!

struct PoleTerm
    term ::OffsetVector{Dict{Variable, Coefficient}}
end

struct Cohesive
    solution ::Dict{RationalComplex, PoleTerm}
end


struct Pole
    taylor_rows::UnitRange{Int64}
    laurent_rows::IntOffsetVector
    coords::RationalComplex
    r_outer::Float64
end

function Pole(upper_term::Int64,
              first_row::Int64,
              coords::RationalComplex,
              r_outer::Float64)

    last_taylor_row = first_row + 2*(upper_term+1)-1
    taylor_rows = first_row : last_taylor_row

    laurent_rows = OffsetArray{Int64}(undef, -upper_term : -1)
    laur_row = last_taylor_row + 1
    for i in eachindex(laurent_rows)
        laurent_rows[i] = laur_row
        laur_row += 2
    end

    return Pole(taylor_rows, laurent_rows, coords, r_outer)
end

function last_row(pole::Pole)
    return maximum(pole.laurent_rows)+1
end

struct EllipticTerm
    re_column::Int64
    # im_column = re_column + 1
    derivative_n::Int64
end

struct SinglePoleElliptic
    terms::Vector{EllipticTerm}
    pole::Pole
    upper_power::Int64
end

function SinglePoleElliptic(upper_term::Int64,
                            first_column::Int64,
                            pole::Pole)

    terms = Vector{EllipticTerm}(undef, 0)

    re_column = first_column
    for i in 1 : upper_term
        push!(terms, EllipticTerm(re_column, i-2))
        re_column += 2
    end

    return SinglePoleElliptic(terms, pole, upper_term)
end

function last_column(elliptic::SinglePoleElliptic)
    return elliptic.terms[end].re_column + 1
end

struct GeneralElliptic
    sp_terms::Vector{SinglePoleElliptic}
end

function fill_slae_coupling_taylor!(slae::Array{Float64, 2},
                                    rows::UnitRange{Int64},
                                    term::EllipticTerm,
                                    delta::RationalComplex,
                                    series_dict::RComplex2OffsetMatrix,
                                    r_source::Float64,
                                    r_dest::Float64,
                                    temp_vector::ComplexOffsetVector)
    upper_term = div(size(rows)[1], 2) - 1

    term_expansion_normalized!(temp_vector, series_dict, delta,
                               term.derivative_n, upper_term, r_source, r_dest)

    first_row = first(rows)
    re_col = term.re_column
    im_col = re_col + 1

    for i in 0 : upper_term
        c, d = reim(temp_vector[i])

        slae[first_row + 2*i,   re_col] =  c
        slae[first_row + 2*i+1, re_col] =  d
    end
    for i in 0 : upper_term
        c, d = reim(temp_vector[i])

        slae[first_row + 2*i,   im_col] = -d
        slae[first_row + 2*i+1, im_col] =  c
    end
end

function fill_slae_coupling_laurent!(slae::Array{Float64, 2},
                                     rows::OffsetArray{Int64, 1, Array{Int64, 1}},
                                     term::EllipticTerm)
    row = rows[term.derivative_n]
    val = (-1.0)^term.derivative_n
    slae[row,   term.re_column]   = val
    slae[row+1, term.re_column]   = 0.0
    slae[row,   term.re_column+1] = 0.0
    slae[row+1, term.re_column+1] = val
end

function doubly_periodic_distance(z1::RationalComplex, z2::RationalComplex)
    delta = z1 - z2

    if real(delta) > 1//2
        delta -= 1
    end
    if real(delta) < -1//2
        delta += 1
    end
    if imag(delta) > 1//2
        delta -= 1im
    end
    if imag(delta) < -1//2
        delta += 1im
    end

    return delta
end

function rational_complex_distance(z1::RationalComplex, z2::RationalComplex,
                                   wei::Weierstrass)
        a, b = reim(z1-z2)
        return abs(a*wei.ω1 + b*wei.ω3)
end

function fill_slae_coupling!(slae::Array{Float64, 2},
                             elliptic::GeneralElliptic,
                             series_dict::RComplex2OffsetMatrix)

    max_n_terms = 0
    for el in elliptic.sp_terms
        # В системе уравнений на каждое комплексное равенство отводится
        # две строки для действительной и мнимой частей.
        n_terms = el.upper_power
        if n_terms > max_n_terms
            max_n_terms = n_terms
        end
    end
    temp_vector = OffsetVector{ComplexF64}(undef, 0 : max_n_terms-1)

    for spe in elliptic.sp_terms
        for term in spe.terms
            fill_slae_coupling_laurent!(slae, spe.pole.laurent_rows, term)
        end
    end

    for source in elliptic.sp_terms
        r_source = source.pole.r_outer
        for term in source.terms
            for dest in elliptic.sp_terms
                delta = doubly_periodic_distance(dest.pole.coords, source.pole.coords)
                r_dest= dest.pole.r_outer
                fill_slae_coupling_taylor!(slae, dest.pole.taylor_rows, term,
                                           delta, series_dict,
                                           r_source, r_dest, temp_vector)
            end
        end
    end
end

function push_to_wei_deltas!(deltas_dict::RComplex2IntDict, delta::RationalComplex,
                              new_value::Int64)
    if real(delta) < 0
        delta = -delta
    end

    old_value = get(deltas_dict, delta, 0)
    if new_value > old_value
        deltas_dict[delta] = new_value
    end
end


#####################################################################
##  TESTING SECTION
#####################################################################



# function test_construct_dict_of_wei_series()
#     ω1 = complex(1.0)
#     ω3 = exp(1im)
#     wei = Weierstrass(ω1, ω3)

#     test_dict = Dict{RationalComplex, Tuple{Int, Bool}}()
#     n = 1000
#     for i in 1 : n
#         a, b, c, d = rand(Int16, 4)
#         r = abs(a) < abs(b) ? a//b : b//a
#         i = abs(c) < abs(d) ? c//d : d//c
#         z = complex(r, i)

#         test_dict[z] = (15, rand(Bool))
#     end

#     return construct_dict_of_wei_series(test_dict, wei)
# end

# function test_term_expansion_normalized!()
#     ω1 = complex(1.0)
#     ω3 = exp(1im)
#     wei = Weierstrass(ω1, ω3)

#     test_dict = Dict{RationalComplex, Tuple{Int, Bool}}()
#     n = 1000
#     n_terms = 15
#     for i in 1 : n
#         a, b, c, d = rand(Int16, 4)
#         r = abs(a) < abs(b) ? a//b : b//a
#         i = abs(c) < abs(d) ? c//d : c//d
#         z = complex(r, i)
#         if r < 0
#             z = -z
#         end

#         test_dict[z] = (n_terms, rand(Bool))
#     end
#     test_dict[complex(0//1)] = (n_terms, rand(Bool))

#     wei_series_dict = construct_dict_of_wei_series(test_dict, wei)

#     output = OffsetVector{ComplexF64}(undef, 0:n_terms)
#     delta = complex(0//1)
#     r_s = 1.0
#     r_d = 1.0
#     @time begin
#         for n_d in -1:13
#             term_expansion_normalized!(output, wei_series_dict, delta, n_d, 12, r_s, r_d)
#         end
#     end
# end

# function test_fill_slae_coupling()

#     ω1 = complex(1.0)
#     ω3 = exp(1im)
#     wei = Weierstrass(ω1, ω3)

#     poles = Vector{Pole}(undef, 0)
#     elliptics = Vector{SinglePoleElliptic}(undef, 0)

#     upper_term = 15
#     first_row = 1
#     first_column = 1
#     central_coords = complex(1//2, 1//2)
#     r_central = 0.4
#     central_fiber = Pole(upper_term, first_row, central_coords, r_central)
#     elliptic = SinglePoleElliptic(upper_term, first_column, central_fiber)
#     push!(elliptics, elliptic)


#     upper_term = 6
#     nx = 20
#     ny = 20
#     r_subfiber = min(abs(ω1)/nx, abs(ω3)/ny) / 2.0

#     for i in 1 : nx
#         for j in 1 : ny
#             coords = complex(i//nx, j//ny)
#             dist = rational_complex_distance(coords, central_coords, wei)

#             if dist > r_central + r_subfiber
#                 first_row = last_row(elliptics[end].pole) + 1
#                 subfiber_pole = Pole(upper_term, first_row, coords, r_subfiber)

#                 first_column = last_column(elliptics[end]) + 1
#                 elliptic = SinglePoleElliptic(upper_term, first_column, subfiber_pole)

#                 push!(elliptics, elliptic)
#             end
#         end
#     end

#     general_elliptic = GeneralElliptic(elliptics)

#     wei_deltas = Dict{RationalComplex, Tuple{Int64, Bool}}()
#     n_poles = size(elliptics)[1]
#     for i in 1 : n_poles
#         for j in i : n_poles
#             delta = elliptics[i].pole.coords - elliptics[j].pole.coords
#             upper_power = max(elliptics[i].upper_power, elliptics[j].upper_power)
#             sigma = false
#             push_to_wei_deltas!(wei_deltas, delta, (upper_power, sigma))
#         end
#     end
# end

# function run_tests()
#     #test_theta()
#     #test_weierstrass()
#     #test_construct_dict_of_wei_series()
#     #test_term_expansion_normalized!()
#     test_fill_slae_coupling()

#     println("Tests passed")
# end

# function bench_weierstrass()
#     ω1 = complex(1.0)
#     ω3 = exp(1im)
#     wei = Weierstrass(ω1, ω3)

#     z_rand = rand(ComplexF64, 10000)
#     output = OffsetVector{ComplexF64}(undef, -2:10)

#     z = z_rand[1]
#     @show @allocated complute_normalized!(z, abs(z), wei, 10, output)

#     println("Benchmark: complute_normalized!()")
#     @time begin
#         for z in z_rand
#             complute_normalized!(z, abs(z), wei, 10, output)
#         end
#     end
# end

# function run_benchmarks()
#     bench_weierstrass()
# end

# run_tests()
# # run_benchmarks()
=#