module Input

using ..Types: RationalComplex

struct LayerData
    E ::Float64 # модуль Юнга
    ν ::Float64 # коэффициент Пуассона
    r ::Float64 # внешний радиус слоя
end

function check(ld::LayerData)
    ld.E > 0         || error("Negative Young's modulus")
    -1. < ld.ν < 0.5 || error("Poisson's ratio out of bounds")
    ld.r > 0         || error("Nonnegative radius")
end

struct FiberData
    layers :: Vector{LayerData}
    r0     :: Float64 # радиус полости
end

function check(fd::FiberData)
    check.(fd.layers)
    fd.r0 >= 0             || error("Negative inner radius")
    fd.r0 < fd.layers[1].r || error("Inner radius is greater than outer radius")
end

struct InclusionData
    fiber     ::FiberData
    center    ::RationalComplex
    max_power ::Int # максимальная степень при z в решениях
end

function check(id::InclusionData)
    check(id.fiber)
    -1 < real(id.center) < 1 || error("Inclusion center out of cell")
    -1 < imag(id.center) < 1 || error("Inclusion center out of cell")
    id.max_power > 2         || error("Solution max power is too low")
end

struct CohesiveData
    E  :: Float64
    ν  :: Float64
end

function check(cd::CohesiveData)
    cd.E > 0         || error("Negative Young's modulus")
    -1. < cd.ν < 0.5 || error("Poisson's ratio out of bounds")    
end

struct CellData
    l1         :: Float64
    l2         :: Float64
    γ          :: Float64
    inclusions :: Vector{InclusionData}
    cohesive   :: CohesiveData
end

function check(cd::CellData)
    cd.l1 > 0      || error("Nonpositive length of cell")
    cd.l2 > 0      || error("Nonpositive length of cell")
    0 < cd.γ < π/2 || error("Invalid cell angle")
    check.(cd.inclusions)
    check(cd.cohesive)
end

end # module Input