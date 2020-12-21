const RationalComplex = Complex{Rational{Int}}

# OffsetArray{ComplexF64}
const ComplexOffsetMatrix = OffsetArray{ComplexF64, 2, Array{ComplexF64, 2}}

# OffsetVector{ComplexF64}
const ComplexOffsetVector = OffsetArray{ComplexF64, 1, Array{ComplexF64, 1}}

# # OffsetVector{Float64}
# const FloatOffsetVector = OffsetArray{Float64, 1, Array{Float64, 1}}

# OffsetVector{Int}
const IntOffsetVector = OffsetArray{Int, 1, Array{Int, 1}}

# ###############################################################################

# function differentiate!(series::ComplexOffsetVector)
#     for n in firstindex(series)+1 : lastindex(series)
#         series[n-1] = n * series[n]
#     end
#     series[end] = 0.0im
# end

# ###############################################################################

# mutable struct CashedVector{N <: Number}
#     array ::OffsetVector{N}
#     last_cashed ::Int
# end

# function CashedVector{N}(range ::UnitRange{Int}) where {N <: Number}
#     array = OffsetVector{N}(undef, range)
#     last_cashed = firstindex(array) - 1
#     CashedVector(array, last_cashed)
# end

# function Base.getindex(v::CashedVector{N} , index::Int) where {N}
#     if index <= v.last_cashed
#        return v.array[index] 
#     else
#         error("Access to undefined values")
#     end
# end

# function Base.setindex!(v::CashedVector{N}, val::N, index::Int) where {N}
#     if index == v.last_cashed+1
#         v.array[index] = val
#         v.last_cashed += 1
#     elseif index <= v.last_cashed
#         v.array[index] = val
#     else
#         error("Access to undefined values")
#     end
# end

# last_cashed(v::CashedVector{N}) where N = v.last_cashed

# not_cashed(v::CashedVector{N}, ind::Int) where N = (ind > v.last_cashed)

# ###############################################################################

# mutable struct BoundedVector{T <: Number}
#     vector ::OffsetArray{T, 1, Array{T, 1}}
#     bottom ::Int
#     top    ::Int
#     function BoundedVector{T}(indices::UnitRange{Int}) where T <: Number
#         n_ind = size(indices, 1)
#         vector = OffsetVector(zeros(T, n_ind), indices)
#         new{T}(vector, first(indices), last(indices))
#     end
# end


# function set_bounds!(bv::BoundedVector{T}, bottom::Int, top::Int) where T
#     if top < bottom
#         error("Верхняя граница меньше нижней")
#     end
#     if bottom < firstindex(bv.vector)
#         error("Нижняя граница выходит за границы массива")
#     end
#     if top > lastindex(bv.vector)
#         error("Верхняя граница выходит за границы массива")
#     end

#     bv.bottom = bottom
#     bv.top    = top

#     return
# end

# function Base.getindex(bv::BoundedVector{T}, key::Int) where T
#     if bv.bottom <= key <= bv.top
#         @inbounds return bv.vector[key]
#     else
#         error("Index out of bounds")
#     end
# end

# function Base.setindex!(bv::BoundedVector{T}, val::T, key::Int) where T
#     if bv.bottom <= key <= bv.top
#         @inbounds bv.vector[key] = val
#     else
#         error("Index out of bounds")
#     end
# end

# Base.firstindex(bv::BoundedVector{T}) where T = bv.bottom
# Base.lastindex(bv::BoundedVector{T}) where T  = bv.top
# Base.eachindex(bv::BoundedVector{T}) where T  = bv.bottom : bv.top
# Base.ndims(bv::BoundedVector{T}) where T = 1

# function Base.fill!(bv::BoundedVector{T}, val::T) where T 
#     for i in eachindex(bv)
#         bv.vector[i] = val
#     end
# end

# ###############################################################################

# struct VarLinForm{T}
#     terms::OffsetArray{T,1,Array{T, 1}}
# end

# Base.getindex(form::VarLinForm, ind::Int) = getindex(form.terms, ind)
# Base.setindex!(form::VarLinForm{T}, val::T, ind::Int) where T = setindex!(form.terms, val, ind)
# Base.similar(form::VarLinForm) = VarLinForm(similar(form.terms))
# Base.firstindex(form::VarLinForm) = firstindex(form.terms)
# Base.first(form::VarLinForm) = first(form.terms)
# Base.lastindex(form::VarLinForm) = lastindex(form.terms)
# Base.eachindex(form::VarLinForm) = firstindex(form.terms) : lastindex(form.terms)
# Base.axes(form::VarLinForm, d) = axes(form.terms, d)
# Base.ndims(form::VarLinForm) = 1 + ndims(first(form.terms))

# ###############################################################################

# "
# Решетка периодов
# "
# struct Lattice
#     ω1::ComplexF64
#     ω3::ComplexF64
# end


# "Преобразование рационального комплексного числа в обыкновенное"
# raw_complex(l::Lattice, rz::RationalComplex) = 2*real(rz)*l.ω1 + 2*imag(rz)*l.ω3