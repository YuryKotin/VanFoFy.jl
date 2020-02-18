module VanFoFy

using OffsetArrays
#using Markdown


include("typesynonims.jl")

include("input.jl")

include("ellipticals.jl")

include("lshear-main.jl")

function construct_problem(cell::CellData)
    wei = Weierstrass(cell)
    #Q = Q_special(wei)
    #lshear = construct_lshear(cell, wei)
    #plane  = construct_plane(cell, wei, Q)
    #lext   = construct_lext(cell, plane)
end

export theta

end # module
