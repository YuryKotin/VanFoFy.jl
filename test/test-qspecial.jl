module TestQSpecial

using Test, OffsetArrays
using VanFoFy.Testing: my_isapprox
using VanFoFy.Types: raw_complex, Lattice, BoundedVector
using VanFoFy.SpecialWeierstrass: Weierstrass
using VanFoFy.SpecialQ: QSpecial
using VanFoFy.FunctionalTerms: EllipticPraecursor
using VanFoFy.FunctionalTerms: add_term_series!, QSpecialTerm

function test()
    ω1 = complex(1.0)
    ω3 = exp(1im)
    lattice = Lattice(ω1, ω3)
    el_praecursor = EllipticPraecursor(lattice, 10)
    wei = el_praecursor.℘
    Q = el_praecursor.Q

    rz = complex(3//20, 7//20)
    z = raw_complex(lattice, rz)

    @test Q.α ≈ 0.93336333346865274496 atol=1e-14
    @test Q.β ≈ -0.015588254984082231674 + 0.024277268726100292157im atol=1e-14
    @test Q.γ1 ≈ 0.15238040677859618355 + 0.44963371772124743897im atol=1e-14
    @test Q.γ3 ≈ -0.21430777320854296764 + 0.42363008252983652691im atol=1e-14
    
    ###########################################################################

    @testset "Q Special function derivatives values" begin    
        ref_array = OffsetVector([
            -0.51663997136492256779 + 0.59933186734118426564im,
            0.03059293447909888428 + 1.19592812759260613831im,
            0.61792409039690798789 + 0.53104522549332222603im,
            0.21949109770876260028 - 0.23037200859738624081im,
            -0.25667024508338398547 - 0.26175052836614570717im,
            -0.34680442894981433621 + 0.03239450713556108613im,
            -0.09690359568209910845 + 0.21090154612274475410im,
            0.12370320077810871562 + 0.08401792981864268650im,
            0.100108873864963227041 - 0.090530908012789471084im,
            -0.027628709900863662741 - 0.107911340307848935272im,
            -0.076637122750858210907 - 0.012691979234158738341im],
            0:10
        )     
        Q_array = OffsetVector(
            [
                Q[rz, n] for n in 0:10
            ],
            0:10
        )   
        @test my_isapprox(Q_array, ref_array, atol=1e-14)
        
        ###################################
        
        ref_array = OffsetVector([
            0.51663997136492256779 - 0.59933186734118426564im, 
            0.03059293447909888428 + 1.19592812759260613831im, 
            -0.61792409039690798789 - 0.53104522549332222603im, 
            0.21949109770876260028 - 0.23037200859738624081im, 
            0.25667024508338398547 + 0.26175052836614570717im, 
            -0.34680442894981433621 + 0.03239450713556108613im, 
            0.09690359568209910845 - 0.21090154612274475410im, 
            0.12370320077810871562 + 0.08401792981864268650im, 
            -0.100108873864963227041 + 0.090530908012789471084im, 
            -0.027628709900863662741 - 0.107911340307848935272im, 
            0.076637122750858210907 + 0.012691979234158738341im],
            0:10
        )
        Q_array = OffsetVector(
            [
                Q[-rz, n] for n in 0:10
            ],
            0:10
        )   
        @test my_isapprox(Q_array, ref_array, atol=1e-14)
    end
            
    ###########################################################################

    @testset "Q special function Taylor series" begin
        ref_array = OffsetVector([
            0.0im, 
            0.0im, 
            0.0im, 
            1.4763729658582473991 + 0.2104518625714126290im, 
            0.0im, 
            0.0018376550954631914927 - 0.0021276761791455488758im, 
            0.0im, 
            -0.014988091444290608037 - 0.050667468042420345242im, 
            0.0im, 
            0.064881835380165342464 + 0.018881015804034253364im, 
            0.0im],
        0:10)

        Q_array = OffsetVector(
            [
                Q.derivs_series[0,n] for n in 0:10
            ],
            0:10
        )
        @test my_isapprox(Q_array, ref_array, atol=1e-14)

        ###############################

        ref_array = OffsetVector([
            0.36909324146456173876 + 0.05261296564285314337im, 
            0.0im, 
            0.0045941377386579740696 - 0.0053191904478638750084im, 
            0.0im, 
            -0.13114580013754276178 - 0.44334034537117805730im, 
            0.0im, 
            1.3625185429834714146 + 0.3965013318847190882im, 
            0.0im, 
            -0.026873960409698735563 + 0.023419258586845129871im, 
            0.0im, 
            -0.02683361029041383614 - 0.18246080727664615306im],
        0:10)

        Q_array = OffsetVector(
            [
                Q.derivs_series[3,n] for n in 0:10
            ],
            0:10
        )
        @test my_isapprox(Q_array, ref_array, atol=1e-14)

    end

    ###########################################################################

    @testset "QSpecial add_term_series" begin
        r = 0.37 * abs(z)
        R = 0.44 * abs(z)
        derivative = 3

        term = QSpecialTerm(derivative, 0//1im, 1.0+0.0im, r)

        coeffs = BoundedVector{ComplexF64}(0:10)
        fill!(coeffs, 0.0im)

        add_term_series!(coeffs, term, point=rz, norm_r=R, praecursor=el_praecursor)

        ref_array = OffsetVector([
            0.00056315410593177144764 - 0.00059107154635268876423im, 
            -0.0014488003776728687610 - 0.0014774765350371930721im, 
            -0.0025839945392610984970 + 0.0002413672449161728543im, 
            -0.0007412698058994620937 + 0.0016133038929863319742im, 
            0.00083272197843991006480 + 0.00056557612335757473203im, 
            0.00053372414223762962120 - 0.00048265982184853852520im, 
            -0.00010802053015565270797 - 0.00042190316636885624810im, 
            -0.00020717254969746117551 - 0.00003431013070774948276im, 
            -0.000042001131111669440601 + 0.000093368639694710144137im, 
            0.000039664194944864604583 + 0.000036827824548616699810im, 
            0.000023949518606323101497 - 0.000010405235374447323730im],
        0:10)

        @test my_isapprox(coeffs, ref_array, atol=1e-14)
    end

    @testset "QSpecial add_term_series at pole" begin
        r = 0.37 * abs(z)
        R = 0.44 * abs(z)
        derivative = 3

        term = QSpecialTerm(derivative, 0//1im, 1.0+0.0im, r)

        coeffs = BoundedVector{ComplexF64}(0:10)
        fill!(coeffs, 0.0im)

        add_term_series!(coeffs, term, point=0//1im, norm_r=R, praecursor=el_praecursor)

        ref_array = OffsetVector([
            0.00049756420364096220176 + 0.00007092605718652539044im, 
            0.0im, 
            9.675126022313664423e-7 - 1.1202066818049628070e-6im, 
            0.0im, 
            -4.314660106551368770e-6 - 0.000014585773237050365825im, 
            0.0im,
            7.0028400502485037052e-6 + 2.0378698118994099433e-6im, 
            0.0im, 
            -2.1577595219354940800e-8 + 1.8803751825948733033e-8im, 
            0.0im, 
            -3.365813594520620680e-9 - 2.2886561254798310949e-8im],
        0:10)

        @test my_isapprox(coeffs, ref_array, atol=1e-14)
    end
end

end # module TestQSpecial
using .TestQSpecial
TestQSpecial.test()