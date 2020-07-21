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
            n_B_vars += incl.max_power * 2
        end
        n_C_vars = 2
        for incl in inclusions
            n_C_vars += (incl.max_power+2) * 2
        end
    
        first_B = first_index
        last_B = first_B + n_B_vars - 1
        B_inds = first_B : last_B
        first_C = last_B + 1
        last_C = first_C + n_C_vars - 1
        C_inds = first_C : last_C
        all_inds = first_B : last_C
        
        ϕ_terms = OffsetVector{EllipticalTerm}(undef, all_inds)
        ψ_terms = OffsetVector{EllipticalTerm}(undef, all_inds)
        
        ψ_terms[first_B]   = ConstTerm(0.0im)
        ψ_terms[first_B+1] = ConstTerm(0.0im)
        
        for c in C_inds
            ϕ_terms[c] = ConstTerm(0.0im)
        end

        B_var = first_B
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
    
        Φ_terms = OffsetVector(
            [ differentiate(ϕ_terms[b]) for b in all_inds ],
            all_inds
        )
        Φ = VarLinForm(Φ_terms)
        
        
        B_var = first_B + 2
        C_var = first_C
        for k in 1 : K
            N_k = inclusions[k].max_power
            ζ_k = inclusions[k].center
            r_k = inclusions[k].radius
            
            for n in -1 : N_k-2
                factor = -(n+2)/(r_k^2)
                ψ_terms[B_var] = QSpecialTerm(n+1, ζ_k, complex(factor,0.0), r_k)
                B_var += 1
                ψ_terms[B_var] = QSpecialTerm(n+1, ζ_k, complex(0.0,factor), r_k)
                B_var += 1
            end
            for n in -1 : N_k    
                ψ_terms[C_var] = WeierstrassTerm(n, ζ_k, 1.0+0.0im, r_k)
                C_var += 1
                ψ_terms[C_var] = WeierstrassTerm(n, ζ_k, 0.0+1.0im, r_k)
                C_var += 1
            end
        end
        
        ψ_terms[C_var] = ZTerm(1.0+0.0im)
        C_var += 1
        ψ_terms[C_var] = ZTerm(0.0+1.0im)
        C_var += 1

        ψ = EllipticalSolution(ψ_terms)
    
        new(E, ν, praecursor, ϕ, Φ, ψ)
    end
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
