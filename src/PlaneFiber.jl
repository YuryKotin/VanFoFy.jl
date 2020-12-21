struct PlaneLayer
    E ::Float64
    ν ::Float64
    r ::Float64
    # polynomial form of solution
    ϕ ::VarPolyForm
    χ ::VarPolyForm
end

function PlaneLayer(E, ν, r, r0, first_index, max_power)
    pow_range_ϕ = -max_power     : max_power + 2
    pow_range_χ = -max_power - 2 : max_power
    
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
    χ = VarPolyForm(var_range, pow_range_χ)
    
    for i in axes(ar,1)
        ϕ[ar[i],i] = 1.0 + 0.0im
        ϕ[ai[i],i] = 0.0 + 1.0im
    end
    for i in axes(br,1)
        χ[br[i],i] = 1.0 + 0.0im
        χ[bi[i],i] = 0.0 + 1.0im
    end
    
    if r0 > 0
        # TODO Hollow layer
    end

    PlaneLayer(E, ν, r, ϕ, χ)    
end

function displacements_forces!(layer::PlaneLayer, buffer::VarPolyFormBox)
    on_circle!(layer.ϕ, layer.r, buffer.f)
    on_circle_conj!(layer.χ, layer.r, buffer.by)
    
    G = layer.E / (2(1+layer.ν))
    κ = 3 - 4*layer.ν

    displacements_forces!(G, κ, buffer)
end

function PlaneLayer(E, ν, r, prev_layer::PlaneLayer, buffer::VarPolyFormBox)
    G = E / (2(1+ν))
    κ = 3 - 4ν
    r1 = prev_layer.r
    
    displacements_forces!(prev_layer, buffer)
    plane_coupling!(G, κ, buffer)
    ϕ = from_circle(buffer.f, r1)
    χ = from_circle_conj(buffer.by, r1)

    PlaneLayer(E, ν, r, ϕ, χ)    
end

struct PlaneFiber
    layers ::Vector{PlaneLayer}
    r0     ::Float64
    buffer ::VarPolyFormBox
end

function PlaneFiber(data::FiberData, first_index, max_power)
    r0 = data.r0
    
    n_layers = size(data.layers, 1)
    layers = Vector{PlaneLayer}(undef, n_layers)
    
    ld = data.layers[1]
    layers[1] = PlaneLayer(
        ld.E, ld.ν, ld.r,
        r0, first_index, max_power)
    
    buffer = VarPolyFormBox(layers[1].ϕ)
    
    if n_layers > 1
        for n in 2 : n_layers
            ld = data.layers[n]
            layers[n] = PlaneLayer(ld.E, ld.ν, ld.r, layers[n-1], buffer)
        end
    end

    displacements_forces!(layers[end], buffer)

    PlaneFiber(layers, r0, buffer)
end

struct PlaneInclusion
    fiber  ::PlaneFiber
    center ::ComplexF64
end

function PlaneInclusion(data::InclusionData, first_index, center)
    fiber = PlaneFiber(data.fiber, first_index, data.max_power)
    PlaneInclusion(fiber, center)
end