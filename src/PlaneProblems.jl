function displacements_forces!(G, κ, buffer::VarPolyFormBox)
    f = buffer.f
    zbF = buffer.zbF
    by = buffer.by
    displ = buffer.displ
    force = buffer.force
    
    @views for p in powers(ϕ)
        @. zbF[:,2-p] = conj(f[:,p]) * p
    end
    
    @views for p in powers(f)
        @. displ[:,p] = (κ*f[:,p] - zbF[:,p] - by[:,p]) / (2G)
        @. force[:,p] =    f[:,p] + zbF[:,p] + by[:,p]
    end
    
    return
end

function plane_coupling!(G, κ, buffer::VarPolyFormBox)
    displ = buffer.displ
    force = buffer.force

    f = buffer.f
    @views for p in powers(ϕ)
        @. f[:,p] = (1 / (1+κ))*(2G*displ[:,p] + force[:,p])
    end
    
    zbF = buffer.zbf
    @views for p in powers(ϕ)
        @. zbF[:,2-p] = conj(f[:,p]) * p
    end
    
    bc = buffer.by
    @views for p in powers(ϕ)
        @. bc[:,p] = force[:,p] - f[:,p] - zbF[:,p]
    end

    return
end


include("PlaneFiber.jl")
include("PlaneCohesive.jl")

struct PlaneProblem
    inclusions               :: Vector{PlaneInclusion}
    cohesive                 :: PlaneCohesive
end
