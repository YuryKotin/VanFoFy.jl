struct FiberData
    E        :: Vector{Float64}
    μ        :: Vector{Float64}
    G        :: Vector{Float64}
    radii    :: Vector{Float64}
    r_inner  :: Float64
    center   :: RationalComplex
    n_terms  :: Int64
end

struct CohesiveData
    E  :: Float64
    μ  :: Float64
    G  :: Float64
end

struct CellData
    l1       :: Float64
    l2       :: Float64
    γ        :: Float64
    fibers   :: Vector{FiberData}
    cohesive :: CohesiveData
end