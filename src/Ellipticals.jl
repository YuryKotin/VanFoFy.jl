module Ellipticals

using OffsetArrays

using ..TypeSynonyms: RationalComplex, ComplexOffsetMatrix, ComplexOffsetVector
using ..Input: FiberData, CellData

mutable struct CashedVector{N <: Number}
    array ::OffsetVector{N}
    last_cashed ::Int
end

function CashedVector{N}(range ::UnitRange{Int}) where {N <: Number}
    array = OffsetVector{N}(undef, range)
    last_cashed = firstindex(array) - 1
    CashedVector(array, last_cashed)
end

function Base.getindex(v::CashedVector{N} , index::Int) where {N}
    if index <= v.last_cashed
       return v.array[index] 
    else
        error("Access to undefined values")
    end
end

function Base.setindex!(v::CashedVector{N}, val::N, index::Int) where {N}
    if index == v.last_cashed+1
        v.array[index] = val
        v.last_cashed += 1
    elseif index <= v.last_cashed
        v.array[index] = val
    else
        error("Access to undefined values")
    end
end

last_cashed(v::CashedVector{N}) where N = v.last_cashed

not_cashed(v::CashedVector{N}, ind::Int) where N = (ind > v.last_cashed)

abstract type EllipticFunction end

include("Ellipticals-Theta.jl")
include("Ellipticals-Weierstrass.jl")
include("Ellipticals-QSpecial.jl")
include("Ellipticals-Series.jl")

struct EllipticPraecursor
    ℘ ::Weierstrass
    Q ::QSpecial
end

function EllipticPraecursor(
                            ω1::ComplexF64,
                            ω3::ComplexF64,
                            max_derivative ::Int)

    ℘ = Weierstrass(ω1, ω3, max_derivative)
    Q = QSpecial(℘)
    EllipticPraecursor(℘, Q)
end
#=
struct Elliptical
    ℘ ::Weierstrass
    Q ::QSpecial
    ℘_series ::RComplex2OffsetMatrix
    Q_series ::RComplex2OffsetMatrix
end

function max_power_per_delta(fibers::Vector{FiberData})
    delta_to_nterms = RComplex2IntDict()
    delta_to_nterms[0] = 0

    n_fibers = length(fibers)
    for i in 1 : n_fibers
        z1 = fibers[i].center
        n1 = fibers[i].n_terms
        
        for j in i+1 : n_fibers
            z2= fibers[j].center
            delta = z2 - z1
            if real(delta) < 0
                delta = -delta
            end

            n_terms = get(delta_to_nterms, delta, -1)
            n2 = fibers[j].n_terms
            n = max(n1, n2)
            if n_terms < n
                delta_to_nterms[delta] = n
            end
        end

        if delta_to_nterms[0] < n1
            delta_to_nterms[0] = n1
        end
    end

    return delta_to_nterms
end

function Elliptical(data::CellData)
    l1 = data.l1
    l2 = data.l2
    γ  = data.γ
    
    ω1 = complex(l1 / 2.0)
    ω3 = (l2 / 2.0) * exp(1.0im * γ)
    tolerance = eps(imag(ω3))
    θ = Theta(ω1, ω3, tolerance)
    ℘ = Weierstrass(ω1, ω3, θ)
    Q = QSpecial(℘)

    Δ_dict = max_power_per_delta(data.fibers)
    max_power = findmax(Δ_dict)[1]
    temp_container = OffsetArray{ComplexF64}(undef, -1 : 2*max_power)

    ℘_series = Dict{RationalComplex, ComplexOffsetMatrix}()
    Q_series = Dict{RationalComplex, ComplexOffsetMatrix}()

    for Δ in keys(Δ_dict)
        upper_term = Δ_dict[Δ]

        # Для вычисления разложений в ряд производных нужно разложить исходную
        # функцию с запасом членов.
        # Всего нужно n производных по n членов в каждой.
        diff_upper_term = 2 * upper_term

        if Δ != 0
            z = real(Δ)*ω1 + imag(Δ)*ω3
            norm_factor = abs(z)

            weierstrass_normalized!(z=z,
                                norm_factor=norm_factor,
                                w=℘,
                                upper_term=diff_upper_term,
                                output=temp_container)
        else
            expand_at_pole!(output=temp_container,
                            weierstrass_data=℘,
                            upper_term=diff_upper_term)
        end
        
        ℘_series_array = OffsetArray{ComplexF64}(undef, 0:upper_term, -1:upper_term-2)
        derivative_series_from_initial!(initial_series=temp_container,
                                        derivatives_series=℘_series_array)
        ℘_series[Δ] = ℘_series_array

        qspecial_normalized!(z=z,
                             norm_factor=norm_factor,
                             Q=Q,
                             upper_term=diff_upper_term,
                             ℘_derivs=temp_container)
        Q_series_array = OffsetArray{ComplexF64}(undef, 0:upper_term, 0:upper_term-2)
        derivative_series_from_initial!(initial_series=temp_container,
                                        derivatives_series=℘_series_array)
        ℘_series[Δ] = Q_series_array
    end
    # TODO
end
=#
end # module Ellipticals