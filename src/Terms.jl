#==============================================================================
Функциональные слагаемые

# Export

- abstract type FunctionalTerm:
        абстрактный тип-родитель всех функциональных слагаемых

## Elliptical

- abstract type EllipticalTerm <: FunctionalTerm:
        абстрактный тип-родитель слагаемых, связанных с эллиптическими функциями

- struct WeierstrassTerm <: EllipticalTerm:
        слагаемое ``a r^{n+2}/(n+1)! ℘^{(n)}(z-ζ_k)``

- struct QSpecialTerm <: EllipticalTerm:
        слагаемое ``a r^{n+3}/(n+1)! Q^{(n)}(z-ζ_k)``

- struct ZTerm <: EllipticalTerm:
        слагаемое ``a z``

- struct ConstTerm <: EllipticalTerm:
        слагаемое ``a``

- differentiate(term <: EllipticalTerm)
        дифференциирование слагаемого

- struct EllipticalPraecursor:
        экземпляр эллиптических функций

- add_term_series!(output, term <: EllipticalTerm; 
                   point, factor, norm_r, power_shift, conjugated, praecursor):
        разложение слагаемого в ряд term и прибавление этого ряда к контейнеру output

## Polynomial

- struct PolynomialTerm <: FunctionalTerm:
        тип для слагаемых в виде полинома

- differentiate(term :: PolynomialTerm)
        дифференциирование слагаемого

- add_term_series!(output, term :: PolynomialTerm; factor, power_shift, conjugated):
        разложение слагаемого в ряд term и прибавление этого ряда к контейнеру output
        
# Description

В модуле представлены структуры и функции для работы с функциональными слагаемыми вида
``a_1 z^{k_1} + ... + a_n z^{k_n}``, 
``b r^{n+2}/(n+1)! ℘^{(n)}(z-ζ_k)``,
``c r^{n+3}/(n+1)! Q^{(n)}(z-ζ_k)``.
Задаются операции дифференциирования, комплексного сопряжения, умножения на ``z^k``
и разложения в ряд по степеням ``z``.
==============================================================================#

"
Функциональные слагаемые
"
module FunctionalTerms

using ..Types: Lattice, RationalComplex, raw_complex
using ..Types: BoundedVector, ComplexOffsetVector
using ..SpecialWeierstrass: Weierstrass
using ..SpecialQ: QSpecial

using OffsetArrays

###############################################################################

abstract type FunctionalTerm end

include("TermsElliptical.jl")
include("TermsPolynomial.jl")

end # module Terms