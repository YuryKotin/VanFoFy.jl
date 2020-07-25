#==============================================================================
# Org

* TODO Полые волокна
==============================================================================#

module PlaneProblems

using ..Types: VarLinForm, BoundedVector, RationalComplex, set_bounds!
using ..FunctionalTerms: PolynomialTerm, PolynomialSolution
using ..FunctionalTerms: EllipticalSolution
using ..FunctionalTerms: z_conj_diff, conjugate, reconjugate, add!
using ..FunctionalTerms: EllipticPraecursor, EllipticalTerm
using ..FunctionalTerms: WeierstrassTerm, QSpecialTerm, ZTerm, ConstTerm
using ..FunctionalTerms: differentiate, add_term_series!
using ..Input: LayerData, InclusionData, CellData
using ..Types: raw_complex, FloatOffsetVector, ComplexOffsetVector
import ..FunctionalTerms: eachpower

using OffsetArrays

include("PlaneFiber.jl")
include("PlaneCohesive.jl")

struct PlaneProblem
    cohesive ::PlaneCohesive
    fibers   ::Vector{PlaneFiber}
    matrix   ::Matrix{Float64}
    vector   ::Vector{Float64}
end

function PlaneProblem(cell::CellData, praecursor::EllipticPraecursor)
    fibers_data = cell.fibers
    n_fibers = size(fibers_data, 1)
    cohesive_data = cell.cohesive

    inclusions = [
        InclusionData(f.center, f.radii[end], f.max_power) 
        for f in fibers_data
    ]
    cohesive = PlaneCohesive(
        cohesive_data.E, 
        cohesive_data.ν, 
        inclusions, 
        1, 
        praecursor
    )

    fibers = Vector{PlaneFiber}(undef, n_fibers)
    index = lastindex(cohesive) + 1
    for f in 1 : n_fibers
        data = fibers_data[f]
        max_power = data.max_power
        layers = [
            LayerData(
                data.E[i],
                data.ν[i],
                data.radii[i]
            )
            for i in eachindex(data.E)
        ]
        vars = index : index + 2*(2*max_power+4) - 1
        fibers[f] = PlaneFiber(layers, vars)
        index = last(vars) + 1
    end

    slae_size = index - 1

    matrix = zeros(Float64, slae_size, slae_size)
    vector = zeros(Float64, slae_size)
    
    maximum_power = 0
    for f in fibers_data
        if f.max_power > maximum_power
            maximum_power = f.max_power
        end
    end

    buffer = BoundedVector{ComplexF64}(-maximum_power : maximum_power + 2)

    poles = [fibers_data[k].center for k in 1 : n_fibers]
    radii = [fibers_data[k].radii[end] for k in 1 : n_fibers]
    
    func = [displacements_coupling!, forces_coupling!]
    
    for v in eachindex(cohesive)
        row = 1

        for k in 1 : n_fibers
            n_terms = fibers_data[k].max_power
            set_bounds!(buffer, -n_terms, n_terms+2)

            for i in 1 : 2
                # fill!(buffer, 0.0im)

                func[i](
                    buffer, 
                    cohesive, 
                    poles[k], 
                    radii[k], 
                    v
                )
                
                for n in -n_terms : n_terms+2
                    coeff = buffer[n]
                    matrix[row, v] = real(coeff)
                    row += 1    
                    matrix[row, v] = imag(coeff)
                    row += 1
                end
            end
        end
    end

    row = 1
    row_k = 1
    for k in 1 : n_fibers
        n_terms = fibers_data[k].max_power
        set_bounds!(buffer, -n_terms, n_terms+2)
        
        fiber = fibers[k]
        f_index = eachindex(fiber)
        for v in f_index
            row_k = row

            for i in 1:2

                # fill!(buffer, 0.0im)

                func[i](
                    buffer,
                    fiber,
                    v
                )

                for n in -n_terms : n_terms+2
                    coeff = -buffer[n]
                    matrix[row_k, v] = real(coeff)
                    row_k += 1    
                    matrix[row_k, v] = imag(coeff)
                    row_k += 1
                end
            end
        end

        row = row_k
    end

    for v in eachindex(cohesive)
        row = slae_size - 3

        matrix[row, v] = cohesive.σ22[v]
        row += 1
        matrix[row, v] = cohesive.σ23[v]
        row += 1
        matrix[row, v] = cohesive.σ33[v]
        row += 1
        matrix[row, v] = cohesive.rotation[v]
    end

    PlaneProblem(cohesive, fibers, matrix, vector)
end
end # module PlaneProblems