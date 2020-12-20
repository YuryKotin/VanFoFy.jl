struct PlaneLayer
    E ::Float64
    ν ::Float64
    r ::Float64
    # polynomial form of solution
    ϕ ::VarPolyForm
    ψ ::VarPolyForm
end

function PlaneLayer(E, ν, r, r0, first_index, max_power)
    pow_range_ϕ = -max_power     : max_power + 2
    pow_range_ψ = -max_power - 2 : max_power
    
    ar = OffsetVector{Int}(undef, 0:max_power+2)
    br = OffsetVector{Int}(undef, 0:max_power)
    ai = similar(ar)
    bi = similar(br)

    ind = first_index
    for i in axes(ar,1)
        ar[i] = ind; ind += 1
        ai[i] = ind; ind += 1
    end
    for i in axes(br,1)
        br[i] = ind; ind += 1
        bi[i] = ind; ind += 1
    end 
    
    var_range = first(ar) : last(bi)
    
    ϕ = VarPolyForm(var_range, pow_range_ϕ)
    ψ = VarPolyForm(var_range, pow_range_ψ)
    
    for i in axes(ar,1)
        ϕ[ar[i],i] = 1.0 + 0.0im
        ϕ[ai[i],i] = 0.0 + 1.0im
    end
    for i in axes(br,1)
        ψ[br[i],i] = 1.0 + 0.0im
        ψ[bi[i],i] = 0.0 + 1.0im
    end
    
    if r0 > 0
        # TODO Hollow layer
    end

    PlaneLayer(E, ν, r, ϕ, ψ)    
end

function displacements_forces(layer::PlaneLayer, tmp)
    for i in 1 : 5
        axes(tmp[i]) == axes(layer.ϕ) || 
            error("Temp container has wrong shape")
    end

    E = layer.ϕ
    ν = layer.ν
    r = layer.r
    ϕ = layer.ϕ
    ψ = layer.ψ
    
    # f = ϕ(z)/. z-> rτ
    f = tmp[1]
    @views for p in powers(ϕ)
        factor = r^p
        @. f[:,p] = ϕ[:,p] * factor
    end
    
    # zbF = z \bar Φ /.{z -> rτ, \bar z -> r/τ}
    zbF = tmp[2]
    @views for p in powers(ϕ)
        factor = p * r^p
        @. zbF[:,2-p] = conj(ϕ[:,p]) * factor
    end
    
    # by = \bar ψ(z)/. \bar z -> r/τ
    by = tmp[3]
    @views for p in powers(ψ)
        factor = r^p
        @. by[:,-p] = ψ[:,p] * factor
    end

    G = E / (2(1+ν))
    κ = 3 - 4ν

    displ = tmp[4]
    force = tmp[5]
    @views for p in powers(f)
        factor = 1 / (2G)
        @. displ[:,p] = factor*(κ*f[:,p] - zbF[:,p] - by[:,p])
        @. force[:,p] =           f[:,p] + zbF[:,p] + by[:,p]
    end
    
    return
end

function PlaneLayer(E, ν, r, prev_layer::PlaneLayer, tmp::Vector{VarPolyForm})

end

struct PlaneFiber
    layers ::Vector{PlaneLayer}
    r0     ::Float64
end

function PlaneFiber(data::FiberData, first_index, max_power)
    r0 = data.r0
    
    n_layers = size(data.layers, 1)
    layers = Vector{PlaneLayer}(undef, n_layers)
    
    ld = data.layers[1]
    layers[1] = PlaneLayer(
        ld.E, ld.ν, ld.r,
        r0, first_index, max_power)
    
    if n_layers > 1
        tmp = [similar(layers[1].ϕ) for _ in 1:5]
        for n in 2 : n_layers
            ld = data.layers[n]
            fill!.(tmp, 0.0im)
            layers[n] = PlaneLayer(ld.E, ld.ν, ld.r, layers[n-1], tmp)
        end
    end

    PlaneFiber(layers, r0)
end

struct PlaneInclusion
    fiber  ::PlaneFiber
    center ::ComplexF64
    displ  ::VarPolyForm
    force  ::VarPolyForm
end