###############################################################################
###  Q special function
###############################################################################

struct QSpecial <: EllipticFunction
    ℘ :: Weierstrass
    α ::ComplexF64
    β ::ComplexF64
    γ1::ComplexF64
    γ3::ComplexF64
end

function QSpecial(℘::Weierstrass)

end

function qspecial_normalized!(;
                            z ::ComplexF64,
                            norm_factor::Float64,
                            Q::QSpecial,
                            upper_term::Int64,
                            ℘_derivs::ComplexOffsetVector)

end