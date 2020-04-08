#=
struct Layer
    solution ::OffsetVector{Dict{Variable, Coefficient}}
end

struct Fiber
    layers ::Vector{Layer}
end
=#