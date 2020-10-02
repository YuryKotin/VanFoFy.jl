"
Слой волокна для плоской задачи
"
struct PlaneLayer
    " модуль Юнга"
    E       ::Float64
    "коэффициент Пуассона"
    ν       ::Float64
    "внутренний радиус слоя"
    r_inner ::Float64
    "внешний радиус слоя"
    r_outer ::Float64
    "функции решения в виде линейных форм от нормированных полиномов"
    ϕ             ::PolynomialForm
    inner_z_bar_Φ ::PolynomialForm
    inner_bar_ψ   ::PolynomialForm
    outer_z_bar_Φ ::PolynomialForm
    outer_bar_ψ   ::PolynomialForm
    """
    Конструктор сплошного ядра
    """
    function PlaneLayer(data::LayerData, r_fiber::Float64, var_indices::UnitRange{Int}) 
    
        E = data.E
        ν = data.ν
        r_outer = data.r
    
        if r_outer <= 0 || r_fiber <= 0
            error("Nonpositive radius")
        end
        if r_outer > r_fiber
            error("Layer is bigger than fiber")
        end
        
        n_vars = size(var_indices, 1)
        if n_vars % 4 != 0
            error("Odd number of variables complex parts")
        end
    
        N = n_vars ÷ 4 - 2
        
        ϕ = VarLinForm(
            OffsetVector(
                [
                    PolynomialTerm(-N : N+2, r_fiber) for v in var_indices
                ], 
                var_indices
            )
        )
        bottom_ϕ = first(var_indices)
        top_ϕ    = bottom_ϕ + 2(N+3) - 1
        
        ψ = VarLinForm(
            OffsetVector(
                [
                    PolynomialTerm(-(N+2) : N, r_fiber) for v in var_indices
                ], 
                var_indices
            )
        )       
        bottom_ψ = top_ϕ + 1
        top_ψ    = bottom_ψ + 2(N+1) - 1

        for i in 0 : N+2
            ϕ[bottom_ϕ + 2i + 0][i] = 1.0 + 0.0im
            ϕ[bottom_ϕ + 2i + 1][i] = 0.0 + 1.0im
        end
        
        for i in 0 : N
            ψ[bottom_ψ + 2i + 0][i] = 1.0 + 0.0im
            ψ[bottom_ψ + 2i + 1][i] = 0.0 + 1.0im
        end
    
        outer_z_bar_Φ = z_conj_diff(ϕ, r_outer)
        outer_bar_ψ = conjugate(ψ, r_outer)
    
        inner_z_bar_Φ = similar(outer_z_bar_Φ)
        inner_bar_ψ =   similar(outer_bar_ψ)
    
    
        new(
            E, ν, 0.0, r_outer, ϕ, 
            inner_z_bar_Φ, inner_bar_ψ, outer_z_bar_Φ, outer_bar_ψ
        )
    
    end
    """
    Конструктор следующего слоя
    """
    function PlaneLayer(data::LayerData, prev_layer::PlaneLayer)
        E = data.E
        ν = data.ν
        r_outer = data.r
    
        r_inner = prev_layer.r_outer
        if r_inner >= r_outer
            error("Layer has nonpositive width")
        end
    
        E1 = prev_layer.E
        ν1 = prev_layer.ν
        G1 = E1/(2*(1+ν1))
        κ1 = 3 - 4*ν1
        
        E2 = E
        ν2 = ν
        κ2 = 3 - 4*ν2
        G2 = E2/(2*(1+ν2))
    
        ϕ1       = prev_layer.ϕ
        z_bar_Φ1 = prev_layer.outer_z_bar_Φ
        bar_ψ1   = prev_layer.outer_bar_ψ
    
        ϕ2 = similar(prev_layer.ϕ)
        
        add!(ϕ2, ϕ1,       (1+κ1*G2/G1)/(1+κ2))
        add!(ϕ2, z_bar_Φ1, (1-G2/G1)   /(1+κ2))
        add!(ϕ2, bar_ψ1,   (1-G2/G1)   /(1+κ2))
        
        z_bar_Φ2 = z_conj_diff(ϕ2, r_inner)
        
        bar_ψ2 = similar(prev_layer.outer_bar_ψ)
    
        add!(bar_ψ2, ϕ1,       +1.0)
        add!(bar_ψ2, z_bar_Φ1, +1.0)
        add!(bar_ψ2, bar_ψ1,   +1.0)
        add!(bar_ψ2, ϕ2,       -1.0)
        add!(bar_ψ2, z_bar_Φ2, -1.0)
    
        outer_z_bar_Φ2 = z_conj_diff(ϕ2, r_outer)
        outer_bar_ψ2 = reconjugate(bar_ψ2, r_inner, r_outer)
    
        new(
            E2, ν2, r_inner, r_outer, ϕ2, 
            z_bar_Φ2, bar_ψ2, outer_z_bar_Φ2, outer_bar_ψ2
        )
    end
    """
    """
    function PlaneLayer(data::LayerData, r_fiber::Float64, ϵ11::Float64) 
    
        E = data.E
        ν = data.ν
        κ = 3 - 4ν
        r_outer = data.r
    
        if r_outer <= 0 || r_fiber <= 0
            error("Nonpositive radius")
        end
        if r_outer > r_fiber
            error("Layer is bigger than fiber")
        end

        ϕ = VarLinForm(
            OffsetVector(
                [
                    PolynomialTerm(1 : 1, r_fiber) for v in 0 : 0
                ], 
                0 : 0
            )
        )
        
        ψ = VarLinForm(
            OffsetVector(
                [
                    PolynomialTerm(-1 : -1, r_fiber) for v in 0 : 0
                ], 
                0 : 0
            )
        )       
        
        ϕ[0][1] = - E * ν / ((κ-1)*(ν+1)) * ϵ11 * r_fiber
        ψ[0][-1] = 0.0
    
        outer_z_bar_Φ = z_conj_diff(ϕ, r_outer)
        outer_bar_ψ = conjugate(ψ, r_outer)
    
        inner_z_bar_Φ = similar(outer_z_bar_Φ)
        inner_bar_ψ =   similar(outer_bar_ψ)
    
    
        new(
            E, ν, 0.0, r_outer, ϕ, 
            inner_z_bar_Φ, inner_bar_ψ, outer_z_bar_Φ, outer_bar_ψ
        )
    
    end
end

Base.firstindex(layer::PlaneLayer) = firstindex(layer.ϕ)
Base.lastindex(layer::PlaneLayer) = lastindex(layer.ϕ)
Base.eachindex(layer::PlaneLayer) = firstindex(layer) : lastindex(layer)
eachpower(layer::PlaneLayer) = eachpower(layer.ϕ)

"
Волокно в плоской задаче
"
struct PlaneFiber
    "Слои волокна"
    layers  :: Vector{PlaneLayer}
    "Конструктор"
    function PlaneFiber(data::Vector{LayerData}, var_indices::UnitRange{Int})
        if data[1].r <= 0
            error("Layer has nonpositive outer radius")
        end
        for n in 2 : lastindex(data)
            if data[n].r <= data[n-1].r
                error("Layer has nonpositive width")
            end
        end
        
        r_fiber = data[end].r
        
        layers = Vector{PlaneLayer}(undef, size(data, 1))
        
        layers[1] = PlaneLayer(data[1], r_fiber, var_indices)
        for n in 2 : lastindex(layers)
            layers[n] = PlaneLayer(data[n], layers[n-1])
        end
        
        new(layers)
    end
end

Base.eachindex(fiber::PlaneFiber) = eachindex(fiber.layers[1])
eachpower(fiber::PlaneFiber) = eachpower(fiber.layers[1])

function displacements_coupling!(output, fiber::PlaneFiber, var::Int)
    layer = fiber.layers[end]
    
    E = layer.E
    ν = layer.ν
    G = E/(2*(1+ν))
    κ = 3 - 4ν

    ϕ = layer.ϕ[var]
    z_bar_Φ = layer.outer_z_bar_Φ[var]
    bar_ψ = layer.outer_bar_ψ[var]

    for p in eachindex(ϕ)
        output[p] = (κ*ϕ[p] - z_bar_Φ[p] - bar_ψ[p]) / (2G)
    end
end

function forces_coupling!(output, fiber::PlaneFiber, var::Int)
    layer = fiber.layers[end]

    E = layer.E
    ν = layer.ν
    G = E/(2*(1+ν))
    κ = 3 - 4ν

    ϕ = layer.ϕ[var]
    z_bar_Φ = layer.outer_z_bar_Φ[var]
    bar_ψ = layer.outer_bar_ψ[var]

    for p in eachindex(ϕ)
        output[p] = ϕ[p] + z_bar_Φ[p] + bar_ψ[p]
    end
end
