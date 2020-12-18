module VanFoFy

include("all_modules.jl")

using ..Input: CellData

struct VanFoFyProblem
    doubly_periodic_functions :: DoublyPeriodicFunction
    longitudinal_shear        :: LongitudinalShear
    plane_problem             :: PlaneProblem
    longitudinal_extension    :: LongitudinalExtension
end

function VanFoFyProblem(cell::CellData)
    doubly_periodic_funcions = DoublyPeriodicFunction(cell)
    longitudinal_shear = LongitudinalShear(cell, doubly_periodic_funcions)
    plane_problem = PlaneProblem(cell, doubly_periodic_funcions)
    longitudinal_extension = LongitudinalExtension(cell, plane_problem)
    VanFoFyProblem(
        doubly_periodic_funcions, 
        longitudinal_shear, 
        plane_problem, 
        longitudinal_extension)
end

struct VanFoFySolution
    longitudinal_shear     :: LongitudinalShearSolution
    plane_problem          :: PlaneProblemSolution
    longitudinal_extension :: LongitudinalExtensionSolution
end

function VanFoFySolution(problem::VanFoFyProblem, σ)
    
end

function compute(solution::VanFoFySolution, z::ComplexF64)
    # TODO проконтролировать, чтобы точка для вычислений принадлежала ячейке периодичности
    
    σ13, σ23 = 
        compute_longitudinal_shear(solution.longitudinal_shear, z)
    σ11_p, σ22_p, σ12_p = 
        compute_plane_problem(solution.plane_problem, z)
    σ33, σ11_e, σ22_e, σ12_e = 
        compute_longitudinal_extension(solution.longitudinal_extension, z)
    
    σ11 = σ11_p + σ11_e
    σ22 = σ22_p + σ22_e
    σ12 = σ12_p + σ12_e
    
    (σ11, σ22, σ12, σ13, σ23, σ33)
end

end # module
