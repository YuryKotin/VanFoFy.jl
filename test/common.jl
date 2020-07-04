module Common

using Test, OffsetArrays
using VanFoFy.Types: VarLinForm, BoundedVector
using VanFoFy.FunctionalTerms: PolynomialTerm

function my_isapprox(
    x::VarLinForm{T}, 
    y::VarLinForm{T}; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : √eps(), 
) where T
    
    x_bottom = firstindex(x)
    y_bottom = firstindex(y)
    if x_bottom != y_bottom
        println("First var indices don't match: $x_bottom, $y_bottom")
        return false
    end
    
    x_top = lastindex(x)
    y_top = lastindex(y)
    if x_top != y_top
        println("Last var indices don't match: $x_top, $y_top")
        return false
    end

    for v in eachindex(x)
        if !my_isapprox(x[v], y[v], atol=atol, rtol=rtol)
            println("Error in test on variable $v")
            return false
        end
    end
    
    return true
end

function my_isapprox(
    x::PolynomialTerm, 
    y::PolynomialTerm; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : √eps(), 
)
    
    x_bottom = firstindex(x)
    y_bottom = firstindex(y)
    if x_bottom != y_bottom
        println("First powers don't match: $x_bottom, $y_bottom")
        return false
    end

    x_top = lastindex(x)
    y_top = lastindex(y)
    if x_top != y_top
        println("Last powers don't match: $x_top, $y_top")
        return false
    end

    for i in eachindex(x)
        if ! isapprox(x[i], y[i], atol=atol, rtol=rtol)
            println("Values don't match: x[$i]=", x[i], " y[$i]=", y[i])
            return false
        end
    end

    return true
end

function my_isapprox(
    x::OffsetArray{T,1,Array{T, 1}}, 
    y::OffsetArray{T,1,Array{T, 1}}; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : √eps(), 
) where T
    
    x_bottom = firstindex(x)
    y_bottom = firstindex(y)
    if x_bottom != y_bottom
        println("First indices don't match: $x_bottom, $y_bottom")
        return false
    end

    x_top = lastindex(x)
    y_top = lastindex(y)
    if x_top != y_top
        println("Last indices don't match: $x_top, $y_top")
        return false
    end

    for i in eachindex(x)
        if ! isapprox(x[i], y[i], atol=atol, rtol=rtol)
            println("Values don't match: x[$i]=", x[i], " y[$i]=", y[i])
            return false
        end
    end

    return true
end

function my_isapprox(
    x::BoundedVector{T}, 
    y::OffsetArray{T,1,Array{T, 1}}; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : √eps(), 
) where T
    
    x_bottom = firstindex(x)
    y_bottom = firstindex(y)
    if x_bottom != y_bottom
        println("First indices don't match: $x_bottom, $y_bottom")
        return false
    end

    x_top = lastindex(x)
    y_top = lastindex(y)
    if x_top != y_top
        println("Last indices don't match: $x_top, $y_top")
        return false
    end

    for i in eachindex(x)
        if ! isapprox(x[i], y[i], atol=atol, rtol=rtol)
            println("Values don't match: x[$i]=", x[i], " y[$i]=", y[i])
            return false
        end
    end

    return true
end


function coincede_indices(x::VarLinForm{T}, y::VarLinForm{T}) where T
    
    x_bottom = firstindex(x)
    y_bottom = firstindex(y)
    if x_bottom != y_bottom
        println("First var indices don't match: $x_bottom, $y_bottom")
        return false
    end
    
    x_top = lastindex(x)
    y_top = lastindex(y)
    if x_top != y_top
        println("Last var indices don't match: $x_top, $y_top")
        return false
    end

    for v in eachindex(x)

        x_v_bottom = firstindex(x[v])
        y_v_bottom = firstindex(y[v])
        if x_v_bottom != y_v_bottom
            println("In var $v first poly indices don't match: $x_v_bottom, $y_v_bottom")
            return false
        end
        
        x_v_top = lastindex(x[v])
        y_v_top = lastindex(y[v])
        if x_v_top != y_v_top
            println("In var $v last poly indices don't match: $x_v_top, $y_v_top")
            return false
        end
    end

    return true
end

end # module Common
