module Ellipticals

using ..Types: Lattice
using ..SpecialWeierstrass: Weierstrass
using ..SpecialQ: QSpecial

struct EllipticPraecursor
    l ::Lattice
    ℘ ::Weierstrass
    Q ::QSpecial
    function EllipticPraecursor(l ::Lattice, max_derivative ::Int)
        ℘ = Weierstrass(l, max_derivative)
        Q = QSpecial(℘)
        new(l, ℘, Q)
    end
end

end # module Ellipticals