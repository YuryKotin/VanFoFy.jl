#==============================================================================
# Org

* TODO Полые волокна
==============================================================================#

module PlaneProblems

using ..FunctionalTerms: PolynomialTerm, PolynomialSolution
using ..FunctionalTerms: z_conj_diff, conjugate, reconjugate, add!
using ..Types: VarLinForm

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
end


"""
Сплошное ядро
"""
function PlaneLayer(E::Float64, 
                    ν::Float64,  
                    r_outer::Float64, 
                    r_fiber::Float64, 
                    var_indices::UnitRange{Int}) 

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


    PlaneLayer(E, ν, 0.0, r_outer, ϕ, 
                inner_z_bar_Φ, inner_bar_ψ, outer_z_bar_Φ, outer_bar_ψ)

end

"""
Следующий слой
"""
function PlaneLayer(E::Float64, ν::Float64, r_outer::Float64, prev_layer::PlaneLayer)
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

    PlaneLayer(E2, ν2, prev_layer.r_outer, r_outer, ϕ2, 
                z_bar_Φ2, bar_ψ2, outer_z_bar_Φ2, outer_bar_ψ2)
end

end # module PlaneProblems