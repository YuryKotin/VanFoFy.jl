module TestSymbolic

using Test, OffsetArrays
#using VanFoFy.Symbolic: VarLinForm, variables, add!, add_conjugated!, mul!
using VanFoFy.Symbolic: EllipticalTerm
using VanFoFy.Symbolic: WeierstrassTerm, QSpecialTerm, ZTerm, ConstTerm

function test()
    @testset "SymbolicSolution" begin
        ss = OffsetVector{EllipticalTerm}(undef, 10:20)

        ss[15] = WeierstrassTerm(0.0im, 1.0+0.0im, 1.0, 1)
        @test ss[15] isa WeierstrassTerm

        ss[14] = QSpecialTerm(0.0im, 1.0+0.0im, 1.0, 1)
        @test ss[14] isa QSpecialTerm
    end
    #=
    @testset "VarLinForm" begin
        form1 = VarLinForm{ComplexF64}()
        @test length(variables(form1)) == 0
        
        @test form1[5] == 0.0im

        form1[5] = 2.0+3.0im
        @test form1[5] == 2.0+3.0im

        form2 = VarLinForm{ComplexF64}()
        form2[3] = 1.0-8.0im
        add!(form1, form2, 2)
        @test form1[5] == 2.0+3.0im
        @test form1[3] == 2.0 - 16.0im

        add_conjugated!(form1, form2, 1.5)
        @test form1[5] == 2.0+3.0im
        @test form1[3] == 3.5 - 4.0im

        add!(form1, form1)
        @test form1[5] == 4.0+6.0im
        @test form1[3] == 7.0 - 8.0im

        add_conjugated!(form1, form1)
        @test form1[5] == 8.0+0.0im
        @test form1[3] == 14.0 - 0.0im

        form1[5] += 4.0im
        mul!(form1, 0.5)
        @test form1[5] == 4.0+2.0im
        @test form1[3] == 7.0 - 0.0im

        mul!(form1, 2)
        @test form1[5] == 8.0+4.0im
        @test form1[3] == 14.0 - 0.0im

        mul!(form1, 1im)
        @test form1[5] == 8.0im-4.0
        @test form1[3] == 14.0im - 0.0
    end
    =#
end

end # module TestSymbolic
using .TestSymbolic
TestSymbolic.test()