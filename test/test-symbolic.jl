module TestSymbolic

using Test
using VanFoFy: VarLinForm, variables, add!, add_conjugated!, mul!

function test()
    @testset "VarLinForm" begin
        form1 = VarLinForm()
        @test length(variables(form1)) == 0
        
        @test form1[5] == 0.0im

        form1[5] = 2.0+3.0im
        @test form1[5] == 2.0+3.0im

        form2 = VarLinForm()
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
end
end # module TestSymbolic
using .TestSymbolic
TestSymbolic.test()