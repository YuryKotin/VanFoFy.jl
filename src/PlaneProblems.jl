#==============================================================================
# Org

* TODO Полые волокна
==============================================================================#

module PlaneProblems

using ..Types: VarLinForm, BoundedVector
using ..FunctionalTerms: PolynomialTerm, PolynomialSolution
using ..FunctionalTerms: EllipticalSolution
using ..FunctionalTerms: z_conj_diff, conjugate, reconjugate, add!
using ..FunctionalTerms: EllipticPraecursor, EllipticalTerm
using ..FunctionalTerms: WeierstrassTerm, QSpecialTerm, ZTerm, ConstTerm
using ..FunctionalTerms: differentiate
using ..Input: LayerData, InclusionData

using OffsetArrays

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
    ϕ             ::PolynomialSolution
    inner_z_bar_Φ ::PolynomialSolution
    inner_bar_ψ   ::PolynomialSolution
    outer_z_bar_Φ ::PolynomialSolution
    outer_bar_ψ   ::PolynomialSolution
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
        if n_vars % 2 == 1
            error("Odd number of variables complex parts")
        end
    
        max_power = n_vars - 1
    
        ϕ = VarLinForm(
            OffsetVector(
                [
                    PolynomialTerm(-max_power:max_power+2, r_fiber) for v in var_indices
                ], 
                var_indices
            )
        )
        
        ψ = VarLinForm(
            OffsetVector(
                [
                    PolynomialTerm(-max_power-2:max_power, r_fiber) for v in var_indices
                ], 
                var_indices)
            )       
    
        bottom_var = first(var_indices)
        for v in var_indices
            p = (v - bottom_var) ÷ 2
            ϕ_poly = ϕ[v]
            ψ_poly = ψ[v]
            if (v - bottom_var) % 2 == 0
                ϕ_poly[p] = 1.0+0.0im
                ψ_poly[p] = 1.0+0.0im
            else
                ϕ_poly[p] = 0.0+1.0im
                ψ_poly[p] = 0.0+1.0im
            end
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
    
        bar_ψ2 = similar(prev_layer.outer_bar_ψ)
        
        add!(ϕ2, ϕ1,       (1+κ1*G2/G1)/(1+κ2))
        add!(ϕ2, z_bar_Φ1, (1-G2/G1)   /(1+κ2))
        add!(ϕ2, bar_ψ1,   (1-G2/G1)   /(1+κ2))
        
        z_bar_Φ2 = z_conj_diff(ϕ2, r_inner)
    
        add!(bar_ψ2, ϕ1,       +1.0)
        add!(bar_ψ2, z_bar_Φ1, +1.0)
        add!(bar_ψ2, bar_ψ1,   +1.0)
        add!(bar_ψ2, ϕ2,       -1.0)
        add!(bar_ψ2, z_bar_Φ2, -1.0)
    
        outer_z_bar_Φ2 = z_conj_diff(ϕ2, r_outer)
        outer_bar_ψ2 = reconjugate(bar_ψ2, r_inner, r_outer)
    
        new(
            E2, ν2, prev_layer.r_outer, r_outer, ϕ2, 
            z_bar_Φ2, bar_ψ2, outer_z_bar_Φ2, outer_bar_ψ2
        )
    end
end


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

Base.getindex(fiber::PlaneFiber, ind::Int) = fiber.layers[ind]
Base.firstindex(fiber::PlaneFiber) = firstindex(fiber.layers)
Base.lastindex(fiber::PlaneFiber) = lastindex(fiber.layers)

function displacements_series!(output, fiber::PlaneFiber, var::Int)
    layer = fiber[end]
    
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

function forces_series!(output, fiber::PlaneFiber, var::Int)
    layer = fiber[end]

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

struct PlaneCohesive
    " модуль Юнга"
    E ::Float64
    "коэффициент Пуассона"
    ν ::Float64
    "прекурсор эллиптических функций"
    praesursor ::EllipticPraecursor
    "функции решения в виде линейных форм эллиптических функций"
    ϕ ::EllipticalSolution
    Φ ::EllipticalSolution
    ψ ::EllipticalSolution
    "Конструктор"
    function PlaneCohesive(
        E ::Float64,
        ν ::Float64,
        inclusions ::Vector{InclusionData},
        first_index ::Int,
        praecursor ::EllipticPraecursor
    )
    
        n_B_vars = 2
        for incl in inclusions
            n_B_vars += (incl.max_power + 2) * 2
        end
    
        B_inds = first_index : first_index + n_B_vars - 1
        ϕ_terms = OffsetVector{EllipticalTerm}(undef, B_inds)
        
        B_var = first_index
        ϕ_terms[B_var] = ZTerm(1.0+0.0im)
        B_var += 1
        ϕ_terms[B_var] = ZTerm(0.0+1.0im)
        B_var += 1
    
        K = size(inclusions, 1)
        for k in 1 : K
            N_k = inclusions[k].max_power
            ζ_k = inclusions[k].center
            r_k = inclusions[k].radius
    
            for n in -1 : N_k-2
                ϕ_terms[B_var] = WeierstrassTerm(n, ζ_k, 1.0+0.0im, r_k)
                B_var += 1
                ϕ_terms[B_var] = WeierstrassTerm(n, ζ_k, 0.0+1.0im, r_k)
                B_var += 1
            end
        end
    
        ϕ = EllipticalSolution(ϕ_terms)
    
        Φ = EllipticalSolution(
            OffsetVector(
                [
                    differentiate(ϕ_terms[b]) for b in B_inds
                ],
                B_inds
            )
        )
        
        C_var = B_var
        n_C_vars = n_B_vars
        
        C_inds = first_index + 2 : C_var + n_C_vars - 1
        ψ_terms = OffsetVector{EllipticalTerm}(undef, C_inds)
    
        ψ_terms[C_var] = ZTerm(1.0+0.0im)
        C_var += 1
        ψ_terms[C_var] = ZTerm(0.0+1.0im)
        C_var += 1

        B_var = first_index + 2
        for k in 1 : K
            N_k = inclusions[k].max_power
            ζ_k = inclusions[k].center
            r_k = inclusions[k].radius
    
            for n in -1 : N_k-2
                ψ_terms[B_var] = QSpecialTerm(n, ζ_k, -1.0-0.0im, r_k)
                B_var += 1
                ψ_terms[B_var] = QSpecialTerm(n, ζ_k, -0.0-1.0im, r_k)
                B_var += 1
    
                ψ_terms[C_var] = WeierstrassTerm(n, ζ_k, 1.0+0.0im, r_k)
                C_var += 1
                ψ_terms[C_var] = WeierstrassTerm(n, ζ_k, 0.0+1.0im, r_k)
                C_var += 1
            end
        end
    
        ψ = EllipticalSolution(ψ_terms)
    
        new(E, ν, praecursor, ϕ, Φ, ψ)
    end
end


end # module PlaneProblems