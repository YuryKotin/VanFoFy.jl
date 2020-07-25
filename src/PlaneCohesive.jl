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
    "осредненные напряжения"
    σ22 ::FloatOffsetVector
    σ23 ::FloatOffsetVector
    σ33 ::FloatOffsetVector
    "кручение ячейки в целом"
    rotation ::FloatOffsetVector
    "Конструктор"
    function PlaneCohesive(
        E ::Float64,
        ν ::Float64,
        inclusions ::Vector{InclusionData},
        first_index ::Int,
        praecursor ::EllipticPraecursor
    )
    
        n_A_vars = 2
        for incl in inclusions
            n_A_vars += incl.max_power * 2
        end
        n_B_vars = 2
        for incl in inclusions
            n_B_vars += (incl.max_power+2) * 2
        end
    
        first_A = first_index
        last_A = first_A + n_A_vars - 1
        A_inds = first_A : last_A
        first_B = last_A + 1
        last_B = first_B + n_B_vars - 1
        B_inds = first_B : last_B
        all_inds = first_A : last_B
        
        ϕ_terms = OffsetVector{EllipticalTerm}(undef, all_inds)
        ψ_terms = OffsetVector{EllipticalTerm}(undef, all_inds)
        
        ψ_terms[first_A]   = ConstTerm(0.0im)
        ψ_terms[first_A+1] = ConstTerm(0.0im)
        
        for b in B_inds
            ϕ_terms[b] = ConstTerm(0.0im)
        end

        A_var = first_A
        ϕ_terms[A_var] = ZTerm(1.0+0.0im)
        A_0_ind = A_var
        A_var += 1
        ϕ_terms[A_var] = ZTerm(0.0+1.0im)
        A_var += 1
    
        K = size(inclusions, 1)
        A_m1_inds = zeros(Int, K)
        for k in 1 : K
            N_k = inclusions[k].max_power
            ζ_k = inclusions[k].center
            r_k = inclusions[k].radius
    
            for n in -1 : N_k-2
                ϕ_terms[A_var] = WeierstrassTerm(n, ζ_k, 1.0+0.0im, r_k)
                if n == -1
                    A_m1_inds[k] = A_var
                end
                A_var += 1
                ϕ_terms[A_var] = WeierstrassTerm(n, ζ_k, 0.0+1.0im, r_k)
                A_var += 1
            end
        end
    
        ϕ = EllipticalSolution(ϕ_terms)
    
        Φ_terms = OffsetVector(
            [ differentiate(ϕ_terms[b]) for b in all_inds ],
            all_inds
        )
        Φ = VarLinForm(Φ_terms)
        
        
        A_var = first_A + 2
        B_var = first_B
        B_m1_inds = zeros(Int, K)
        for k in 1 : K
            N_k = inclusions[k].max_power
            ζ_k = inclusions[k].center
            r_k = inclusions[k].radius
            
            for n in -1 : N_k-2
                factor = -(n+2)/(r_k^2)
                ψ_terms[A_var] = QSpecialTerm(n+1, ζ_k, complex(factor,0.0), r_k)
                A_var += 1
                ψ_terms[A_var] = QSpecialTerm(n+1, ζ_k, complex(0.0,factor), r_k)
                A_var += 1
            end
            for n in -1 : N_k    
                ψ_terms[B_var] = WeierstrassTerm(n, ζ_k, 1.0+0.0im, r_k)
                if n == -1
                    B_m1_inds[k] = B_var
                end
                B_var += 1
                ψ_terms[B_var] = WeierstrassTerm(n, ζ_k, 0.0+1.0im, r_k)
                B_var += 1
            end
        end
        
        ψ_terms[B_var] = ZTerm(1.0+0.0im)
        B_0_ind = B_var
        B_var += 1
        ψ_terms[B_var] = ZTerm(0.0+1.0im)
        B_var += 1

        ψ = EllipticalSolution(ψ_terms)

        ω1 = praecursor.l.ω1
        ω3 = praecursor.l.ω3
        η1 = praecursor.℘.η1
        η3 = praecursor.℘.η3
        γ1 = praecursor.Q.γ1
        γ3 = praecursor.Q.γ3
        e_i_θ = ω3 / abs(ω3)
        e_2i_θ = e_i_θ^2
        cosθ = real(e_i_θ)
        sinθ = imag(e_i_θ)
        l2 = 2 * abs(ω1)
        l3 = 2 * abs(ω3)
        κ = 3 - 4ν

        radii = [incl.radius for incl in inclusions]

        Σ11 = c_Σ11(all_inds, A_0_ind, A_m1_inds, praecursor, radii)

        Σ12 = c_Σ12(all_inds, A_0_ind, A_m1_inds, praecursor, radii)

        Σ21 = c_Σ21(all_inds, A_m1_inds, B_0_ind, B_m1_inds, praecursor, radii)

        Σ22 = c_Σ22(all_inds, A_m1_inds, B_0_ind, B_m1_inds, praecursor, radii)

        rotation = imag(
            shift(all_inds, A_0_ind, A_m1_inds, B_0_ind, B_m1_inds, praecursor, radii, ν)
        )

        ###########################

        σ332 = zero_vector(Float64, all_inds)
        σ222 = zero_vector(Float64, all_inds)
        σ232 = zero_vector(Float64, all_inds)
        σ333 = zero_vector(Float64, all_inds)
        σ223 = zero_vector(Float64, all_inds)
        σ233 = zero_vector(Float64, all_inds)

        # σ332 = (2 Re[Σ11] + Re[Σ21]) / l2;
        @__dot__ σ332 = (2 * real(Σ11) + real(Σ21)) / l2
        
        # σ222 = (2 Re[Σ11] - Re[Σ21]) / l2;
        @__dot__ σ222 = (2 * real(Σ11) - real(Σ21)) / l2

        # σ232 = Im[Σ21] / l2;
        @__dot__ σ232 = imag(Σ21) / l2

        # σ333 = (2 Re[Σ12] + Re[Σ22]) / l3;
        @__dot__ σ333 = (2 * real(Σ12) + real(Σ22)) / l3

        # σ223 = (2 Re[Σ12] - Re[Σ22]) / l3;
        @__dot__ σ223 = (2 * real(Σ12) - real(Σ22)) / l3

        # σ233 = Im[Σ22]/l3;
        @__dot__ σ233 = imag(Σ22) / l3

        σ22 = zero_vector(Float64, all_inds)
        σ23 = zero_vector(Float64, all_inds)
        σ33 = zero_vector(Float64, all_inds)

        # σ22 = σ222 l2 Sin[θ] - σ232 l2 Cos[θ] + σ233 l3 Cos[θ];
        @__dot__ σ22 += σ222 * l2 * sinθ
        @__dot__ σ22 -= σ232 * l2 * cosθ
        @__dot__ σ22 += σ233 * l3 * cosθ

        # σ23 = σ232 l2 Sin[θ] - σ332 l2 Cos[θ] + σ233 l3 Sin[θ] + σ333 l3 Cos[θ];
        @__dot__ σ23 += σ232 * l2 * sinθ
        @__dot__ σ23 -= σ332 * l2 * cosθ
        @__dot__ σ23 += σ233 * l3 * sinθ
        @__dot__ σ23 += σ333 * l3 * cosθ

        # σ33 = σ333 l3 Sin[θ];
        @__dot__ σ33 += σ333 * l3 * sinθ

        new(E, ν, praecursor, ϕ, Φ, ψ, σ22, σ23, σ33, rotation)
    end
end

function zero_vector(T, inds)
    OffsetVector(
        zeros(T, size(inds, 1)),
        inds
    ) 
end

function c_Σ11(inds, A_0, A_m1, praecursor, r)
    Σ11 = zero_vector(ComplexF64, inds)

    ω1 = praecursor.l.ω1
    η1 = praecursor.℘.η1
    
    # A[0] 2 ω1 + Sum[A[k, -1] r[k] (-2 η1), {k, Nfibers}]
    Σ11[A_0 + 0] = (2ω1) * (1.0+0.0im)
    Σ11[A_0 + 1] = (2ω1) * (0.0+1.0im)
    for k in axes(r,1)
        Σ11[A_m1[k] + 0] = r[k] * (- 2η1) * (1.0+0.0im)
        Σ11[A_m1[k] + 1] = r[k] * (- 2η1) * (0.0+1.0im)
    end
    return Σ11
end

function c_Σ12(inds, A_0, A_m1, praecursor, r)
    Σ12 = zero_vector(ComplexF64, inds)

    ω3 = praecursor.l.ω3
    η3 = praecursor.℘.η3
    e_i_θ = ω3 / abs(ω3)
    
    # (A[0] 2 ω3 + Sum[A[k, -1] r[k] (-2 η3), {k, Nfibers}]) Exp[-I θ]
    Σ12[A_0 + 0] = (2ω3) / e_i_θ * (1.0+0.0im)
    Σ12[A_0 + 1] = (2ω3) / e_i_θ * (0.0+1.0im)
    for k in axes(r, 1)
        Σ12[A_m1[k] + 0] = r[k] * (- 2η3) / e_i_θ * (1.0+0.0im)
        Σ12[A_m1[k] + 1] = r[k] * (- 2η3) / e_i_θ * (0.0+1.0im)
    end

    return Σ12
end

function c_Σ21(inds, A_m1, B_0, B_m1, praecursor, r)
    Σ21 = zero_vector(ComplexF64, inds)

    ω1 = praecursor.l.ω1
    η1 = praecursor.℘.η1
    γ1 = praecursor.Q.γ1
    
    # B[0] 2 ω1 +
    # Sum[B[k, -1] r[k] (-2 η1), {k, Nfibers}] -
    # Sum[A[k, -1] r[k] (γ1 - 2 η1), {k, Nfibers}];
    Σ21[B_0 + 0] = (2ω1) * (1.0+0.0im)
    Σ21[B_0 + 1] = (2ω1) * (0.0+1.0im)
    for k in axes(r,1)
        Σ21[B_m1[k] + 0] = r[k] * (- 2η1) * (1.0+0.0im)
        Σ21[B_m1[k] + 1] = r[k] * (- 2η1) * (0.0+1.0im)

        Σ21[A_m1[k] + 0] = -r[k] * (γ1 - 2η1) * (1.0+0.0im)
        Σ21[A_m1[k] + 1] = -r[k] * (γ1 - 2η1) * (0.0+1.0im)
    end    

    return Σ21
end

function c_Σ22(inds, A_m1, B_0, B_m1, praecursor, r)
    Σ22 = zero_vector(ComplexF64, inds)

    ω3 = praecursor.l.ω3
    η3 = praecursor.℘.η3
    γ3 = praecursor.Q.γ3
    e_i_θ = ω3 / abs(ω3)
    e_2i_θ = e_i_θ^2

    # (B[0] 2 ω3 + Sum[B[k, -1] r[k] (-2 η3), {k, Nfibers}] -
    # Sum[A[k, -1] r[k] (γ3 - 2 η3 Exp[-2 I θ]), {k, Nfibers}]  ) Exp[-I θ]

    Σ22[B_0 + 0] = (2ω3) / e_i_θ * (1.0+0.0im)
    Σ22[B_0 + 1] = (2ω3) / e_i_θ * (0.0+1.0im)
    for k in axes(r, 1)
        Σ22[B_m1[k] + 0] = r[k] * (-2η3) / e_i_θ * (1.0+0.0im)
        Σ22[B_m1[k] + 1] = r[k] * (-2η3) / e_i_θ * (0.0+1.0im)

        Σ22[A_m1[k] + 0] = -r[k] * (γ3 - 2η3 / e_2i_θ) / e_i_θ * (1.0+0.0im)
        Σ22[A_m1[k] + 1] = -r[k] * (γ3 - 2η3 / e_2i_θ) / e_i_θ * (0.0+1.0im)
    end

    return Σ22
end

function shift(inds, A_0, A_m1, B_0, B_m1, praecursor, r, ν)
    s = zero_vector(ComplexF64, inds)

    ω1 = praecursor.l.ω1
    η1 = praecursor.℘.η1
    γ1 = praecursor.Q.γ1
    κ = 3 - 4ν
    
    #=
    κ (A[0] 2 ω1 + Sum[A[k, -1] r[k] (-2 η1), {k, Nfibers}]) -
    Conjugate[2 ω1 A[0]] +
    Conjugate[
        -2 ω1 B[0] +
        Sum[B[k, -1] r[k] (2 η1), {k, Nfibers}] +
        Sum[A[k, -1] r[k] γ1, {k, Nfibers}]
    ]
    =#

    s[A_0 + 0] = κ * (2ω1) * (1.0 + 0.0im) - conj((2ω1) * (1.0 + 0.0im))
    s[A_0 + 1] = κ * (2ω1) * (0.0 + 1.0im) - conj((2ω1) * (0.0 + 1.0im))
    s[B_0 + 0] = -conj((2ω1) * (1.0 + 0.0im))
    s[B_0 + 1] = -conj((2ω1) * (0.0 + 1.0im))
    for k in axes(r,1)
        s[A_m1[k] + 0] = κ * r[k] * (-2η1) * (1.0+0.0im)
        s[A_m1[k] + 0] += conj(r[k] * γ1 * (1.0+0.0im))
        s[A_m1[k] + 1] = κ * r[k] * (-2η1) * (0.0+1.0im)
        s[A_m1[k] + 1] += conj(r[k] * γ1 * (0.0+1.0im))

        s[B_m1[k] + 0] = conj(r[k] * (2η1) * (1.0+0.0im))
        s[B_m1[k] + 1] = conj(r[k] * (2η1) * (0.0+1.0im))
    end

    return s
end

Base.firstindex(coh::PlaneCohesive) = firstindex(coh.ϕ)
Base.lastindex(coh::PlaneCohesive) = lastindex(coh.ψ)
Base.eachindex(coh::PlaneCohesive) = firstindex(coh) : lastindex(coh)

function forces_coupling!(
    output, 
    cohesive::PlaneCohesive,
    pole::RationalComplex,
    norm_r::Float64,
    var::Int
)
    R = norm_r

    lattice = cohesive.praesursor.l
    c = raw_complex(lattice, pole)

    ϕ = cohesive.ϕ[var]
    Φ = cohesive.Φ[var]
    ψ = cohesive.ψ[var]

    fill!(output, 0.0im)

    add!(f, factor; power_shift=0, conjugated=false) = 
    add_term_series!(
        output, 
        f, 
        point=pole, 
        factor=complex(factor),
        norm_r=R,
        power_shift=power_shift,
        conjugated=conjugated,
        praecursor=cohesive.praesursor
    )
    
    # ϕ 
    add!(ϕ, 1.0)

    # z bar Φ
    add!(Φ, c, conjugated=true)
    add!(Φ, R, power_shift=1, conjugated=true)

    # bar ψ
    add!(ψ, 1.0, conjugated=true)
end

function displacements_coupling!(
    output, 
    cohesive::PlaneCohesive,
    pole::RationalComplex,
    norm_r::Float64,
    var::Int
)
    E = cohesive.E
    ν = cohesive.ν
    G = E/(2*(1+ν))
    κ = 3 - 4ν

    R = norm_r

    lattice = cohesive.praesursor.l
    c = raw_complex(lattice, pole)

    ϕ = cohesive.ϕ[var]
    Φ = cohesive.Φ[var]
    ψ = cohesive.ψ[var]

    fill!(output, 0.0im)
    
    add!(f, factor; power_shift=0, conjugated=false) = 
    add_term_series!(
        output, 
        f, 
        point=pole, 
        factor=complex(factor),
        norm_r=R,
        power_shift=power_shift,
        conjugated=conjugated,
        praecursor=cohesive.praesursor
    )

    # ϕ 
    add!(ϕ, κ/(2G))

    # z bar Φ
    add!(Φ, -c/(2G), conjugated=true)
    add!(Φ, -R/(2G), power_shift=1, conjugated=true)

    # bar ψ
    add!(ψ, -1/(2G), conjugated=true)
    
end
