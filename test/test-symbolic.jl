module TestSymbolic

using Test, OffsetArrays
#using VanFoFy.Symbolic: VarLinForm, variables, add!, add_conjugated!, mul!
using VanFoFy.Testing: my_isapprox
using VanFoFy.FunctionalTerms: EllipticalTerm, differentiate
using VanFoFy.FunctionalTerms: WeierstrassTerm, QSpecialTerm, ZTerm, ConstTerm
using VanFoFy.FunctionalTerms: add_term_series!, fill!
using VanFoFy.FunctionalTerms: EllipticPraecursor
using VanFoFy.Types: raw_complex, Lattice, differentiate!, BoundedVector
using VanFoFy.Types: set_bounds!
using VanFoFy.FunctionalTerms: PolynomialTerm, conjugate, z_conj_diff

function test()
    @testset "Ellipticals complex computing check" begin
        ω1 = complex(1.0)
        ω3 = exp(1im)
        lattice = Lattice(ω1, ω3)
        el_praecursor = EllipticPraecursor(lattice, 10)
        wei = el_praecursor.℘
        Q = el_praecursor.Q

        rz = complex(3//10, 7//10)
        z = raw_complex(lattice, rz)

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

        #######################################################################
        ## Прямое вычисление разложения в ряд вне полюса
        #######################################################################

        set_bounds!(coeffs, 0, 10)

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
        
        @test my_isapprox(coeffs, ref_array, atol=1e-14)

        set_bounds!(coeffs, -10, 10)


        #######################################################################
        ## Прямое вычисление разложения в ряд около полюса
        #######################################################################

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

        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ## Разложения в ряд производной вне полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        diff_terms = Vector{EllipticalTerm}(undef, 10)
        for n in 1:10
            diff_terms[n] = differentiate(terms[n])
        end

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=rz, norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            1.11564950673748297660 - 0.39489123408373782986im,
            0.21330721451221834606 + 0.59767301830264984996im, 
            0.10956126993447068418 - 0.55921335004997774210im, 
            -0.30052575863632435826 + 0.33161003313314041385im, 
            -0.13592193692920667702 - 0.28891870085433229987im, 
            0.25069638972673385924 + 0.36077206496732300289im, 
            -0.24446082219715611905 - 0.09294747019717171377im, 
            0.21160806112174254667 - 0.09860380017150485732im, 
            -0.08820514147929334192 + 0.12541356092819092027im, 
            -0.018447385973376157625 - 0.097915262662077659495im, 
            0.048385277643970804606 + 0.046418708590132078118im], 
            -10:10)
        
        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ## Разложения в ряд производной около полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=0//1im, norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            2.6595946760554767252 + 3.7234325464776669712im, 
            -17.711462923785660450im,
            6.7700186465435550076 + 0.0im, 
            7.1562959867367306188im, 
            2.1275474555163254031 - 0.0im, 
            0.0im, 
            0.47795745174522819010 - 0.12466510020323501307im, 
            0.040227988964436411923 - 0.029070909337064181066im, 
            0.007478554564484870409 + 0.013861400957836167436im, 
            0.003611894055632630489 - 0.027215117022347534020im, 
            0.0039238802673843761748 + 0.0004509269432378789559im, 
            -0.00013320291957249284884 - 0.00024360072468681545155im, 
            0.000028599122157742358442 - 0.000024906695921860620029im, 
            0.000069025797963986673320 + 0.000027734681116819736760im, 
            0.000016057994195909370806 + 0.000012895680976961388818im,
            2.1947084954494754229e-6 - 9.7774208409021146283e-6im,
            4.8082636354407545414e-7 - 2.7620103491109986283e-7im], 
            -10:10)
        
        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ## Комплексное сопряжение разложения в ряд производной вне полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=rz, conjugated=true, 
                            norm_r=R, praecursor=el_praecursor)
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

        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ## Комплексное сопряжение разложения в ряд производной около полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=0//1im, conjugated=true, 
                            norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            4.8082636354407545414e-7 + 2.7620103491109986283e-7im,
            2.1947084954494754229e-6 + 9.7774208409021146283e-6im, 
            0.000016057994195909370806 - 0.000012895680976961388818im, 
            0.000069025797963986673320 - 0.000027734681116819736760im, 
            0.000028599122157742358442 + 0.000024906695921860620029im, 
            -0.00013320291957249284884 + 0.00024360072468681545155im, 
            0.0039238802673843761748 - 0.0004509269432378789559im, 
            0.003611894055632630489 + 0.027215117022347534020im, 
            0.007478554564484870409 - 0.013861400957836167436im, 
            0.040227988964436411923 + 0.029070909337064181066im, 
            0.47795745174522819010 + 0.12466510020323501307im, 
            0.0im, 
            2.1275474555163254031 + 0.0im, 
            -7.1562959867367306188im, 
            6.7700186465435550076 - 0.0im, 
            17.711462923785660450im, 
            2.6595946760554767252 - 3.7234325464776669712im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im], -10:10)

        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ##Умноженное на z комплексное сопряжение разложения в ряд производной около полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=0//1im, 
                            conjugated=true, factor=R+0.0im, power_shift=1,
                            norm_r=R, praecursor=el_praecursor)
            #=add_term_series!(coeffs, dt, point=0//1im, 
                            conjugated=true, factor=0.0im, power_shift=0,
                            norm_r=R, praecursor=el_praecursor)=#
        end

        ref_array = OffsetVector([
            0.0im, 
            1.9004570694985840482e-7 + 1.0916793445572139488e-7im,
            8.674543643826681014e-7 + 3.8645070169601703457e-6im,
            6.3468917067368283945e-6 - 5.0969934131779686120e-6im, 
            0.000027282315543501732162 - 0.000010962080034514864568im, 
            0.000011303748714093312725 + 9.844324257443225375e-6im, 
            -0.000052648200966665680315 + 0.000096282723757845617820im, 
            0.0015509062229972885028 - 0.0001782280178622649459im, 
            0.0014275942653626408439 + 0.0107567233129908339240im, 
            0.0029558844874784735886 - 0.0054786924014118519588im, 
            0.015900036232552833965 + 0.011490221700647588673im, 
            0.18891177501033251662 + 0.04927364407697869264im, 
            0.0im, 
            0.84090909090909082835 + 0.0im, 
            -2.8285123966942147256im, 
            2.6758370116453793486 - 0.0im, 
            7.0004220809712434104im, 
            1.0511997443016882769 - 1.4716796420223634545im, 
            0.0im, 
            0.0im, 
            0.0im], -10:10)

        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ##Умноженное на z комплексное сопряжение разложения в ряд производной вне полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=rz, 
                            conjugated=true, factor=R+0.0im, power_shift=1,
                            norm_r=R, praecursor=el_praecursor)
            add_term_series!(coeffs, dt, point=rz, 
                            conjugated=true, factor=z, power_shift=0,
                            norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.060157454751560986150 - 0.002981342217215410262im, 
            -0.051062039077184345759 + 0.037194303899843475514im, 
            0.006759265191915687747 - 0.098311559078877258644im, 
            0.05057156426599195703 + 0.14194820122429627207im, 
            -0.13690730405806933456 - 0.04198366363311439170im, 
            0.28590798482179147211 - 0.06027487688179146430im, 
            -0.16327825696555386958 - 0.02670851129343584929im, 
            -0.06221479245141114250 - 0.28772580231650318749im, 
            -0.37386977545084076890 + 0.31273159380014547981im, 
            0.54001846641830331119 - 0.05867648566339955796im, 
            0.60835306319563919875 + 0.68874137834476689690im, 
            0.44095834856765364629 + 0.15608001024857839911im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im, 
            0.0im], -10:10)
        
        @test my_isapprox(coeffs, ref_array, atol=1e-10)
        
        #######################################################################
        ## Имитация выражения ϕ(z) + z bar Φ(z) вне полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=rz, 
                            conjugated=true, factor=R+0.0im, power_shift=1,
                            norm_r=R, praecursor=el_praecursor)
            add_term_series!(coeffs, dt, point=rz, 
                            conjugated=true, factor=z, power_shift=0,
                            norm_r=R, praecursor=el_praecursor)
        end

        for t in terms
            add_term_series!(coeffs, t, point=rz, 
                            norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0.060157454751560986150 - 0.002981342217215410262im, 
            -0.051062039077184345759 + 0.037194303899843475514im, 
            0.006759265191915687747 - 0.098311559078877258644im, 
            0.05057156426599195703 + 0.14194820122429627207im, 
            -0.13690730405806933456 - 0.04198366363311439170im, 
            0.28590798482179147211 - 0.06027487688179146430im, 
            -0.16327825696555386958 - 0.02670851129343584929im, 
            -0.06221479245141114250 - 0.28772580231650318749im, 
            -0.37386977545084076890 + 0.31273159380014547981im, 
            0.54001846641830331119 - 0.05867648566339955796im, 
            0.8304207656067313925 + 1.6606911329281046630im, 
            0.88191669713530740360 + 0.0im, 
            0.042154635699133899407 + 0.118114562657264582368im, 
            0.014434627950062560095 - 0.073676004828234109567im, 
            -0.029695558826948242892 + 0.032767058940954239821im, 
            -0.010744577482527042575 - 0.022838913552984931032im, 
            0.016514545314661201875 + 0.023765745576396702515im, 
            -0.013803239380744470807 - 0.005248187294123805684im, 
            0.0104547105309888213598 - 0.0048716205922583657312im, 
            -0.0038736571822559321028 + 0.0055077190841060601201im, 
            -0.0007291294268571879610 - 0.0038700821595215768867im], -10:10)

        @test my_isapprox(coeffs, ref_array, atol=1e-10)

        #######################################################################
        ## Имитация выражения ϕ(z) + z bar Φ(z) около полюса
        #######################################################################

        fill!(coeffs, 0.0im)

        for dt in diff_terms
            add_term_series!(coeffs, dt, point=0//1im, 
                            conjugated=true, factor=R+0.0im, power_shift=1,
                            norm_r=R, praecursor=el_praecursor)
            #=add_term_series!(coeffs, dt, point=0//1im, 
                            conjugated=true, factor=0.0im, power_shift=0,
                            norm_r=R, praecursor=el_praecursor)=#
        end

        for t in terms
            add_term_series!(coeffs, t, point=0//1im, 
                            norm_r=R, praecursor=el_praecursor)
        end

        ref_array = OffsetVector([
            0, 
            1.9004570694985840482e-7 + 1.0916793445572139488e-7im,
            8.674543643826681014e-7 + 3.8645070169601703457e-6im,
            6.3468917067368283945e-6 - 5.0969934131779686120e-6im, 
            0.000027282315543501732162 - 0.000010962080034514864568im, 
            -0.21022864511162361878 - 0.29432608408021532220im, 
            -0.0000526482009666657 + 1.7502018029665684651im, 
            -0.89039476432546238449 - 0.00017822801786226481im, 
            0.00142759426536264084 - 1.40349947503411653926im, 
            -0.83795320642161241764 - 0.00547869240141185196im, 
            0.31750541258497860797 + 0.48075550908536851535im, 
            0.37782355002066503324 + 0.0im, 
            0.0079500181162764083087 - 0.0057451108503237873978im, 
            0.84189438573825037260 + 0.00182623080047062216im, 
            0.0003568985663406602 - 2.8312015775224623759im, 
            2.6761471928899789319 + 0.0000356456035724525im, 
            -8.7747001611110e-6 + 7.0004060338506173267im, 
            1.0512013591229332210 - 1.4716810483544002786im,
            3.4102894429377148261e-6 + 1.3702600043143597650e-6im,
            7.0521018963742516541e-7 + 5.6633260146421892289e-7im,
            8.674543643826683132e-8 - 3.8645070169601710868e-7im], -10:10)

    @test my_isapprox(coeffs, ref_array, atol=1e-10)

    end

    @testset "Polynomials, differentiate" begin
        top = 3
        bottom = -3
    
        poly_coeffs = OffsetVector([1.0+0.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_poly = differentiate(poly)
        d_bottom = -4
        d_top = 2

        @test firstindex(d_poly) == d_bottom
        @test lastindex(d_poly)  ==  d_top

        ref_coeffs = OffsetVector([(1.0+0.0im)*(n+1) for n in d_bottom:d_top], d_bottom:d_top)

        @test my_isapprox(d_poly, ref_coeffs, atol=1e-10)

        ################

        top = 3
        bottom = 0
    
        poly_coeffs = OffsetVector([1.0+0.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_poly = differentiate(poly)
        d_bottom = 0
        d_top = 2

        ref_coeffs = OffsetVector([(1.0+0.0im)*(n+1) for n in d_bottom:d_top], d_bottom:d_top)

        @test my_isapprox(d_poly, ref_coeffs, atol=1e-10)

        ################

        top = 0
        bottom = -3
    
        poly_coeffs = OffsetVector([1.0+0.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_poly = differentiate(poly)
        d_bottom = -4
        d_top = -2

        ref_coeffs = OffsetVector([(1.0+0.0im)*(n+1) for n in d_bottom:d_top], d_bottom:d_top)

        @test my_isapprox(d_poly, ref_coeffs, atol=1e-10)

        ################

        top = 0
        bottom = 0
    
        poly_coeffs = OffsetVector([1.0+0.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_poly = differentiate(poly)
        d_bottom = 0
        d_top = 0

        ref_coeffs = OffsetVector([0.0im for n in d_bottom:d_top], d_bottom:d_top)

        @test my_isapprox(d_poly, ref_coeffs, atol=1e-10)
        
    end

    @testset "Polynomials, conjugate" begin
        top = 3
        bottom = -5
    
        poly_coeffs = OffsetVector([1.0+1.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        c_poly = conjugate(poly, 1.0)
        c_bottom = -3
        c_top = 5

        ref_coeffs = OffsetVector([(1.0-1.0im) for n in c_bottom:c_top], c_bottom:c_top)

        @test my_isapprox(c_poly, ref_coeffs, atol=1e-10)
        
    end

    @testset "Polynomials, z_conj_diff" begin
        top = 3
        bottom = -3
    
        poly_coeffs = OffsetVector([1.0+1.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_bottom = -4
        d_top = 2

        c_bottom = -2
        c_top = 4

        z_bottom = -1
        z_top = 5

        z_poly = z_conj_diff(poly, 1.0)

        ref_coeffs = OffsetVector(
            [ 
                3-3im,
                2-2im,
                1-1im,
                0im,
                -1+1im,
                -2+2im,
                -3+3im
            ], 
            z_bottom:z_top
        )

        @test my_isapprox(z_poly, ref_coeffs, atol=1e-10)
        
        ################

        top = 3
        bottom = 0
    
        poly_coeffs = OffsetVector([1.0+1.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_bottom = 0
        d_top = 2

        c_bottom = -2
        c_top = 0

        z_bottom = -1
        z_top = 1

        z_poly = z_conj_diff(poly, 1.0)

        ref_coeffs = OffsetVector(
            [ 3-3im, 2-2im, 1-1im], 
            z_bottom:z_top
        )

        @test my_isapprox(z_poly, ref_coeffs, atol=1e-10)
        
        ################

        bottom = -3
        top = 0
        
        poly_coeffs = OffsetVector([1.0+1.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_bottom = -4
        d_top = -2

        c_bottom = 2
        c_top = 4

        z_bottom = 3
        z_top = 5

        z_poly = z_conj_diff(poly, 1.0)

        ref_coeffs = OffsetVector(
            [ -1+1im, -2+2im, -3+3im], 
            z_bottom:z_top
        )

        @test my_isapprox(z_poly, ref_coeffs, atol=1e-10)
        
        ################

        bottom = 0
        top = 0
        
        poly_coeffs = OffsetVector([1.0+1.0im for n in bottom:top], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, 1.0)

        d_bottom = 0
        d_top = 0

        c_bottom = 0
        c_top = 0

        z_bottom = 0
        z_top = 0

        z_poly = z_conj_diff(poly, 1.0)

        ref_coeffs = OffsetVector([ 0im,], z_bottom:z_top)

        @test my_isapprox(z_poly, ref_coeffs, atol=1e-10)
        
    end

    @testset "Random polynomials" begin
        r = 0.63
        R = 0.47

        top = 3
        bottom = -3
    
        poly_coeffs = OffsetVector([
                 0.20848564666417512825 + 0.20424886648029794145im,
                 0.54304809515214014226 - 0.49182451491269629784im,
                 -0.38469937722687719273 - 0.28893744360020923168im,
                 0.99499588339694211570 + 0.32067493901531696210im,
                 -0.63396242992094098412 - 0.82408093224932121856im,
                 -0.56139952476137144899 - 0.50784056637236263398im,
                 -0.50518022781614257966 - 0.61569524282236232082im 
        ], bottom:top) 
        poly = PolynomialTerm(poly_coeffs, r)

        d_poly = differentiate(poly)
        d_bottom = -4
        d_top = 2

        ref_coeffs = OffsetVector([
                -0.99278879363892913457 - 0.97261364990618059956im,
                -1.7239622068321911463 + 1.5613476663895122787im,
                0.61063393210615413143 + 0.45863086285747484139im,
                0.0im,
                -1.0062895713030808320 - 1.3080649718243193558im,
                -1.7822207135281631363 - 1.6121922741979761717im,
                -2.4056201324578214695 - 2.9318821086779154328im
               
        ], d_bottom:d_top)

        @test my_isapprox(d_poly, ref_coeffs, atol=1e-10)

        ################

        c_poly = conjugate(poly, R)
        c_bottom = -3
        c_top = 3

        ref_coeffs = OffsetVector([
                -0.087094392719408844639 + 0.106147470390227968706im,
                -0.17390067774478032425 + 0.15731010587510854681im,
                -0.35284026396960399552 + 0.45865325758093983266im,
                0.99499588339694211570 - 0.32067493901531696210im,
                -0.69120499240084920523 + 0.51914563768638777308im,
                1.7531095709034989483 + 1.5877456748225193817im,
                1.2092951473640940474 - 1.1847202291441174538im
        ], c_bottom:c_top)

        @test my_isapprox(c_poly, ref_coeffs, atol=1e-10)
        
        #####################

        z_bottom = -1
        z_top = 5

        z_poly = z_conj_diff(poly, R)

        ref_coeffs = OffsetVector([ 
            -0.46945809602082450018 + 0.57215841101695086302im,
            -0.62490881844185897087 + 0.56529091011163945524im,
            -0.63396242992094098412 + 0.82408093224932121856im,
            0.0im,
            1.2419160773376960449 - 0.9327700479752255580im,
            -6.2997663077555356992 - 5.7055342538438953426im,
            -6.5183690899340289704 + 6.3859048295242217819im
        ], z_bottom:z_top)

        @test my_isapprox(z_poly, ref_coeffs, atol=1e-10)
        
    end
end

end # module TestSymbolic
using .TestSymbolic
TestSymbolic.test()