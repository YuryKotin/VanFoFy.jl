module TestWeierstrass

using Test, OffsetArrays
using VanFoFy.Types: Lattice, raw_complex, BoundedVector
using VanFoFy.SpecialWeierstrass: Weierstrass
using VanFoFy.SpecialQ: QSpecial
using VanFoFy.FunctionalTerms: WeierstrassTerm, add_term_series!
using VanFoFy.FunctionalTerms: EllipticPraecursor

function test()
    ω1 = complex(1.0)
    ω3 = exp(1im)
    lattice = Lattice(ω1, ω3)
    
    el_praecursor = EllipticPraecursor(lattice, 20)
    wei = el_praecursor.℘
    Q = el_praecursor.Q
    

    rz = complex(3//10, 7//10)
    z = raw_complex(lattice, rz)


    @testset "Weierstrass construction" begin
        @test wei.g2 ≈ -1.0321084791761130361 - 2.2551981702100225569im atol=1e-13
        @test wei.g3 ≈ 13.8271498100001721809 + 1.9710124059856877210im atol=1e-13
        @test wei.η1 ≈ 0.91777507848457062778 + 0.02427726872609958786im atol=1e-13
    end

    @testset "Weierstrass computation" begin
        abs_z = abs(z)

        ref_array = OffsetVector(
                [0.10922980957434628990 - 0.71842740909634417168im,
                -0.8598534582241565 + 0.7026363517078360im,
                -0.0639302522525613 - 1.2338575450881468im,
                0.315957874556992 + 1.608024389830146im,
                -1.432308324486529 + 0.345694829053258im,
                1.963873759456078 - 0.492648412322271im,
                -0.683501620549705 + 1.453552820019121im,
                -0.898424448464008 - 2.273035157032416im,
                1.667509374595122 + 1.332241142926466im,
                -2.633101865120501 + 0.633539568271906im,
                2.077985541274223 - 2.171576392262631im,
                0.162780226883512 + 3.098964550714613im,
                -2.479278381517072 - 2.856788007534111im],
                -2:10)
        for n in 1:10  ref_array[n] *= (n+1)  end
        for n in -2:10  ref_array[n] *= (abs_z^(n+2))  end

        for n in 10:-1:-1
            @test wei[rz, n] ≈ ref_array[n] atol=1e-10
        end
    end

    @testset "Weierstrass series expansions" begin
        r = 0.37 * abs(z)
        R = 0.44 * abs(z)
    
        coeffs = BoundedVector{ComplexF64}(-10:10)
        fill!(coeffs, 0.0im)

        term = WeierstrassTerm(0, 0, 1.0+0.0im, r)

        add_term_series!(coeffs, term, point=rz, norm_r=R, praecursor=el_praecursor)

        ref_series_0 = OffsetVector(
            [  -0.007062266411672984100 - 0.136302147894561409558im, 
                0.02759094970140508682 +  0.14042036496366341880im, 
                -0.074154076823552683706 + 0.017897459976231100409im, 
                0.053582208760999762431 - 0.013441388453674393486im, 
                -0.009213537293326539215 + 0.019593754736521876159im, 
                -0.005744066363879477694 - 0.014532624097382221273im, 
                0.0049161231398375179746 + 0.0039276909685587836976im, 
                -0.0035065809394849473297 + 0.0008437036956069802291im, 
                0.0012304983696016682897 - 0.0012859190581788635874im, 
                0.00004233186991736170212 + 0.00080590233071269013962im, 
                -0.00028031967035285472543 - 0.00032300280537675245263im],
            0:10
        )
        for n in 0:10
            @test coeffs[n] ≈ ref_series_0[n]
        end

        fill!(coeffs, 0.0im)
        term_4 = WeierstrassTerm(4, 0, 1.0+0.0im, r)
        add_term_series!(coeffs, term_4, point=rz, norm_r=R, praecursor=el_praecursor)

        ref_series_4 = OffsetVector(
            [  -0.0009214092844579307214 + 0.0019594936186668951428im, 
            -0.0028722063577332934865 - 0.0072667501875549505089im, 
            0.0073746293530483368259 + 0.0058918916964707662709im, 
            -0.0122737733187416580422 + 0.0029531409902620719271im, 
            0.0086140079571586352508 - 0.0090019761692170623124im, 
            0.0005334137224169422905 + 0.0101549816478491532201im, 
            -0.0058870680302061141609 - 0.0067834679129252263563im, 
            0.0062260580012766311708 + 0.0008464160569543714577im, 
            -0.0038832757904932869823 + 0.0026856895350562824849im, 
            0.0008983187116395736048 - 0.0030265867901243404675im, 
            0.0009320447700863913930 + 0.0018708927349205316335im],
            0:10
        )

        for n in 0:10
            @test coeffs[n] ≈ ref_series_4[n]
        end

        fill!(coeffs, 0.0im)
        term_m1 = WeierstrassTerm(-1, 0, 1.0+0.0im, r)
        add_term_series!(coeffs, term_m1, point=rz, norm_r=R, praecursor=el_praecursor)

        ref_series_m1 = OffsetVector(
            [ -0.28578753432989584260 + 0.23353364293010461794im, 
            -0.00839837086793544102 - 0.16208904073947844049im, 
            0.016405429552186805847 + 0.083493189978394471984im, 
            -0.029394408830957823725 + 0.007094488639226743075im, 
            0.015929845847864793768 - 0.003996088459200494540im, 
            -0.0021913277886830692318 + 0.0046601362616592572480im, 
            -0.0011384636036517886432 - 0.0028803399111928733138im, 
            0.00083517149865965556368 + 0.00066725251975516025778im, 
            -0.00052124851803154613385 + 0.00012541541421184836482im, 
            0.00016258837316058080040 - 0.00016991122690651651859im, 
            5.034060206388960567*10^-6 + 0.000095837033922590211658im],
            0:10
        )

        for n in 0:10
            @test coeffs[n] ≈ ref_series_m1[n]
        end
    end

    @testset "Weierstrass series expansions at pole" begin

    r = 1.0
    R = 1.0

    coeffs = BoundedVector{ComplexF64}(-10:10)
    fill!(coeffs, 0.0im)

    term = WeierstrassTerm(0, 0, 1.0+0.0im, r)

    add_term_series!(coeffs, term, point=0//1im, norm_r=R, praecursor=el_praecursor)

    ref_series_0 = OffsetVector(
        [   0.0im, 
            0.0im, 
            -0.051605423958805651807 - 0.112759908510501127843im, 
            0.0im, 
            0.49382677892857757789 + 0.07039330021377456146im, 
            0.0im, 
            -0.0033505590617761707625 + 0.0038793485894936988143im, 
            0.0im, 
            -0.004785435871899049176 - 0.016177237776385416715im, 
            0.0im, 
            0.018471565761370376530 + 0.005375340001129344245im],
        0:10
    )
    for n in 0:10
        @test coeffs[n] ≈ ref_series_0[n]
    end

    #######################################################################

    r = 0.37 * abs(z)
    R = 0.44 * abs(z)

    #coeffs = BoundedVector{ComplexF64}(-10:10)
    fill!(coeffs, 0.0im)

    term1 = WeierstrassTerm(0, 0, 1.0+0.0im, r)

    add_term_series!(coeffs, term1, point=0//1im, norm_r=R, praecursor=el_praecursor)

    ref_series_1 = OffsetVector(
        [  0.0im, 
            0.0im, 
            -0.0008905793320548119904 - 0.0019459513419365762933im, 
            0.0im, 
            0.00133134766833253198083 + 0.00018977900775080317332im, 
            0.0im, 
            -1.4111517097587130193e-6 + 1.6338614821826384846e-6im, 
            0.0im, 
            -3.1486001553185400245e-7 - 1.06438900737246364231e-6im, 
            0.0im, 
            1.8986254826586904308e-7 + 5.525117705745383757e-8im],
        0:10
    )
    
    for n in 0:10
        @test coeffs[n] ≈ ref_series_1[n]
    end

    fill!(coeffs, 0.0im)
    term_4 = WeierstrassTerm(4, 0, 1.0+0.0im, r)
    add_term_series!(coeffs, term_4, point=0//1im, norm_r=R, praecursor=el_praecursor)

    ref_series_4 = OffsetVector(
        [  0.000133142794497780546890 + 0.000018979045090907650627im, 
        0.0im, 
        -2.1168551975561177993e-6 + 2.4509399993828030394e-6im, 
        0.0im, 
        -2.2041530051441401216e-6 - 7.4511723099532586561e-6im, 
        0.0im, 
        3.9873539256895752173e-6 + 1.1603446796179551156e-6im, 
        0.0im, 
        -2.7859793973179674329e-8 + 2.4278361257042633198e-8im, 
        0.0im, 
        -2.789038286510059390e-9 - 1.8964655585830648380e-8im],
        0:10
    )

    for n in 0:10
        @test coeffs[n] ≈ ref_series_4[n]
    end

    fill!(coeffs, 0.0im)
    term_m1 = WeierstrassTerm(-1, 0, 1.0+0.0im, r)
    add_term_series!(coeffs, term_m1, point=0//1im, norm_r=R, praecursor=el_praecursor)

    ref_series_m1 = OffsetVector(
        [ 0, 
        0, 
        0, 
        -0.00035302243793163717075 - 0.00077136809950639057248im, 
        0, 
        0.00031664485084665635424 + 0.00004513662887046130333im, 
        0, 
        -2.3973233679298599188e-7 + 2.7756720160631695794e-7im, 
        0, 
        -4.160312517538011280e-8 - 1.4063998896212731118e-7im, 
        0],
        0:10
    )

    for n in 0:10
        @test coeffs[n] ≈ ref_series_m1[n]
    end

    end
end

end # module TestWeierstrass
using .TestWeierstrass
TestWeierstrass.test()