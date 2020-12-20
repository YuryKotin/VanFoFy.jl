module VarPolyForms

using Lazy: @forward
using ..Types: ComplexOffsetMatrix

struct VarPolyForm
    data ::ComplexOffsetMatrix
end

function VarPolyForm(var_range, pow_range)
    data = OffsetArray{ComplexF64, 2}(undef, var_range, pow_range)
    fill!(data, 0.0im)
    VarPolyForm(data)
end

@forward VarPolyForm.data getindex, setindex!, fill!, axes

Base.similar(form::VarPolyForm) = VarPolyForm(similar(form.data))

variables(form::VarPolyForm) = axes(form.data, 1)
powers(form::VarPolyForm) =    axes(form.data, 2)

end #module