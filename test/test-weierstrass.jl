module TestWeierstrass

using Test, OffsetArrays
using VanFoFy.Ellipticals: Weierstrass, raw_complex

function test()
    ω1 = complex(1.0)
    ω3 = exp(1im)
    wei = Weierstrass(ω1, ω3, 10)


    @testset "Weierstrass construction" begin
        @test wei.g2 ≈ -1.0321084791761130361 - 2.2551981702100225569im atol=1e-13
        @test wei.g3 ≈ 13.8271498100001721809 + 1.9710124059856877210im atol=1e-13
        @test wei.η1 ≈ 0.91777507848457062778 + 0.02427726872609958786im atol=1e-13
    end

    @testset "Weierstrass computation" begin
        rz = complex(3//10, 7//10)
        z = raw_complex(wei, rz)
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
        for n in 1:10  ref_array[n] *= n+1  end

        for n in 10:-1:-1
            @test wei[rz, n] ≈ ref_array[n]*(abs_z^(n+2)) atol=1e-10
        end
    end
end

end # module TestWeierstrass
using .TestWeierstrass
TestWeierstrass.test()