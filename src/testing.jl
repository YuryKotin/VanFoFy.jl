module Testing

using OffsetArrays
using ..Types: VarLinForm, BoundedVector, ComplexOffsetMatrix
using ..FunctionalTerms: PolynomialTerm, PolynomialForm

function my_isapprox(
    x::PolynomialForm, 
    y::ComplexOffsetMatrix; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : √eps(), 
)
    y_axes = axes(y)
    x_bottom = firstindex(x)
    y_bottom = first(y_axes[2])
    if x_bottom != y_bottom
        println("First var indices don't match: $x_bottom, $y_bottom")
        return false
    end
    
    x_top = lastindex(x)
    y_top = last(y_axes[2])
    if x_top != y_top
        println("Last var indices don't match: $x_top, $y_top")
        return false
    end
    
    flag = true
    for v in eachindex(x)
        a = firstindex(x[v])
        b = first(y_axes[1])
        if a != b
            println("Bottom powers in var $v don't match: $a, $b")
            return false
        end
        a = lastindex(x[v])
        b = last(y_axes[1])
        if a != b
            println("Top powers in var $v don't match: $a, $b")
            return false
        end

        for p in y_axes[1]
            a = x[v][p]
            b = y[p,v]
            if !isapprox(a, b, atol=atol, rtol=rtol)
                println("Error on variable $v, power $p: $a, $b")
                flag = false
            end
        end
    end
    
    return flag
end

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
    
    flag = true
    for v in eachindex(x)
        if !my_isapprox(x[v], y[v], atol=atol, rtol=rtol)
            println("Error in test on variable $v")
            flag =  false
        end
    end
    
    return flag
end

function my_isapprox(
    x, 
    y; 
    atol::Real=0, 
    rtol::Real=atol>0 ? 0 : √eps(), 
)
    x_dims = ndims(x)
    y_dims = ndims(y)
    if x_dims != y_dims
        error("Not matching dimensions")
    end
    if x_dims == 1
        x_bottom = firstindex(x)
        y_bottom = firstindex(y)
        if x_bottom != y_bottom
            println(typeof(x), " , ", typeof(y))
            println("First powers don't match: $x_bottom, $y_bottom")
            return false
        end

        x_top = lastindex(x)
        y_top = lastindex(y)
        if x_top != y_top
            println(typeof(x), " , ", typeof(y))
            println("Last powers don't match: $x_top, $y_top")
            return false
        end

        flag = true
        for i in eachindex(x)
            if ! isapprox(x[i], y[i], atol=atol, rtol=rtol)
                if flag # Печатает только первый раз
                    println(typeof(x), " , ", typeof(y))
                end
                println("Values don't match: x[$i]=", x[i], "; y[$i]=", y[i])
                flag =  false
            end
        end

        return flag
    end

    if x_dims == 2
        flag = true
        
        x_axes = axes(x)
        y_axes = axes(y)
        for i in 1 : x_dims
            flag = flag && (first(x_axes[i]) == first(y_axes[i]))
            flag = flag && (last(x_axes[i]) == last(y_axes[i]))
        end    
        if !flag
            println("Axes don't match:", x_axes, " , ", y_axes)
            return flag
        end

        flag = true
        for i in eachindex(x)
            if ! isapprox(x[i], y[i], atol=atol, rtol=rtol)
                if flag # Печатает только первый раз
                    println(typeof(x), " , ", typeof(y))
                end
                println("Values don't match: x[$i]=", x[i], "; y[$i]=", y[i])
                flag =  false
            end
        end

        return flag
    end

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

end # module Teting
