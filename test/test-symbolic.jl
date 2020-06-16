module TestSymbolic

using Test, OffsetArrays
#using VanFoFy.Symbolic: VarLinForm, variables, add!, add_conjugated!, mul!
using VanFoFy.Symbolic: EllipticalTerm, differentiate
using VanFoFy.Symbolic: WeierstrassTerm, QSpecialTerm, ZTerm, ConstTerm
using VanFoFy.Symbolic: add_term_series!, BoundedVector, fill!
using VanFoFy.Ellipticals: EllipticPraecursor, raw_complex

function test()
    @testset "SymbolicSolution" begin
        ss = OffsetVector{EllipticalTerm}(undef, 10:20)

        ss[15] = WeierstrassTerm(0.0im, 1.0+0.0im, 1.0, 1)
        @test ss[15] isa WeierstrassTerm

        ss[14] = QSpecialTerm(0.0im, 1.0+0.0im, 1.0, 1)
        @test ss[14] isa QSpecialTerm
    end

    @testset "Complex computing check" begin
        ω1 = complex(1.0)
        ω3 = exp(1im)
        el_praecursor = EllipticPraecursor(ω1, ω3, 10)
        wei = el_praecursor.℘
        Q = el_praecursor.Q

        rz = complex(3//10, 7//10)
        z = raw_complex(wei, rz)

        r = 0.37 * abs(z)
        R = 0.44 * abs(z)

        terms = Vector{EllipticalTerm}(undef, 10)
        terms[1] = ConstTerm(0.3+0.47im)
        terms[2] = ZTerm(0.48-0.12im)
        terms[3] = WeierstrassTerm(-1, 0, 1.0+0.0im, r)
        terms[4] = WeierstrassTerm(0,  0, 0.0-2.0im, r)
        terms[5] = WeierstrassTerm(1,  0, 1.5+0.0im, r)
        terms[6] = WeierstrassTerm(2,  0, 0.0+3.5im, r)
        terms[7] = WeierstrassTerm(3,  0, 0.5+0.7im, r)
        terms[8] = QSpecialTerm(0, 0, -0.5+0.7im, r)
        terms[9] = QSpecialTerm(1, 0, 2.5-3.7im, r)
        terms[10] = QSpecialTerm(4, 0, -1.2-1.1im, r)

        coeffs = BoundedVector{ComplexF64}(-10:10)
        fill!(coeffs, 0.0im)

        for t in terms
            add_term_series!(coeffs, t, point=rz, norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.22206770241109222153 + 0.97194975458333776608im, 
            0.44095834856765381282 - 0.15608001024857837136im, 
            0.042154635699133899407 + 0.118114562657264582368im, 
            0.014434627950062560095 - 0.073676004828234109567im, 
            -0.029695558826948242892 + 0.032767058940954239821im, 
            -0.010744577482527042575 - 0.022838913552984931032im, 
            0.016514545314661201875 + 0.023765745576396702515im, 
            -0.013803239380744470807 - 0.005248187294123805684im, 
            0.0104547105309888213598 - 0.0048716205922583657312im, 
            -0.0038736571822559321028 + 0.0055077190841060601201im, 
            -0.0007291294268571879610 - 0.0038700821595215768867im], 0:10)
        
        for n in 0:10
            @test coeffs[n] ≈ ref_array[n]
        end

        #########################################

        fill!(coeffs, 0.0im)

        for t in terms
            add_term_series!(coeffs, t, point=0//1im, norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            -0.21023994886033769980 - 0.29433592840447275751im, 
            1.7501055202428106305im,
            -0.89194567054845963483+0.0im, 
            -1.4142561983471073628im, 
            -0.84090909090909093937+0.0im, 
            0.30160537635242579135 + 0.46926528738472095270im, 
            0.18891177501033251662 - 0.04927364407697868570im, 
            0.0079500181162764135129 - 0.0057451108503237891326im, 
            0.0009852948291594911239 + 0.0018262308004706182592im, 
            0.0003568985663406601568 - 0.0026891808282477084810im, 
            0.00031018124459945774392 + 0.00003564560357245299868im, 
            -8.774700161110952931e-6 - 0.000016047120626307608617im, 
            1.6148212448704726414e-6 - 1.4063320367776030804e-6im,
            3.4102894429377148261e-6 + 1.3702600043143597650e-6im, 
            7.0521018963742516541e-7 + 5.6633260146421892289e-7im, 
            8.674543643826683132e-8 - 3.8645070169601710868e-7im], -10:10)

        for n in -10:10
            @test coeffs[n] ≈ ref_array[n]
        end

        #########################################

        fill!(coeffs, 0.0im)

        diff_terms = Vector{EllipticalTerm}(undef, 10)
        for n in 1:10
            diff_terms[n] = differentiate(terms[n])
        end

        for t in diff_terms
            add_term_series!(coeffs, t, point=rz, norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.048385277643970804606 - 0.046418708590132078118im, 
            -0.018447385973376157625 + 0.097915262662077659495im, 
            -0.08820514147929334192 - 0.12541356092819092027im, 
            0.21160806112174254667 + 0.09860380017150485732im, 
            -0.24446082219715611905 + 0.09294747019717171377im, 
            0.25069638972673385924 - 0.36077206496732300289im, 
            -0.13592193692920667702 + 0.28891870085433229987im, 
            -0.30052575863632435826 - 0.33161003313314041385im, 
            0.10956126993447068418 + 0.55921335004997774210im, 
            0.21330721451221834606 - 0.59767301830264984996im, 
            1.11564950673748297660 + 0.39489123408373782986im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im], -10:10)

        @test coeffs[-10] ≈ ref_array[-10] atol=1e-14
        
        #=for n in -10:10
            @test coeffs[n] ≈ ref_array[n]
        end=#

    end
end

end # module TestSymbolic
using .TestSymbolic
TestSymbolic.test()