module Input

using ..Types: RationalComplex

struct FiberData
    E        :: Vector{Float64}
    ν        :: Vector{Float64}
    radii    :: Vector{Float64}
    r_inner  :: Float64
    center   :: RationalComplex
    max_power:: Int64
end

"
Данные для конструирования слоя
"
struct LayerData
    "Модуль Юнга"
    E ::Float64
    "Коэффициент Пуассона"
    ν ::Float64
    "Внешний радиус слоя"
    r ::Float64
end

struct CohesiveData
    E  :: Float64
    ν  :: Float64
end

struct InclusionData
    center    ::RationalComplex
    radius    ::Float64
    max_power ::Int
end

struct CellData
    l1       :: Float64
    l2       :: Float64
    γ        :: Float64
    fibers   :: Vector{FiberData}
    cohesive :: CohesiveData
end

end # module Input