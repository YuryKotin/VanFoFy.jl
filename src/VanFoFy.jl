module VanFoFy

include("all_modules.jl")

using .Input: CellData
using .Ellipticals: Weierstrass, theta

function construct_problem(cell::CellData)
    wei = Weierstrass(cell)
    #Q = Q_special(wei)
    #lshear = construct_lshear(cell, wei)
    #plane  = construct_plane(cell, wei, Q)
    #lext   = construct_lext(cell, plane)
end

export theta

end # module
