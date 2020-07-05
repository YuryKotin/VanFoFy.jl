module TestPlaneProblems

using Test

using VanFoFy.Testing: my_isapprox, coincede_indices
using VanFoFy.PlaneProblems: PlaneLayer, add!, PlaneFiber
using VanFoFy.PlaneProblems: forces_series!, displacements_series!
using VanFoFy.Input: LayerData, InclusionData
using VanFoFy.FunctionalTerms: EllipticPraecursor
using VanFoFy.Types: Lattice
using VanFoFy.PlaneProblems: PlaneCohesive

function test1()
    @testset "Fiber layers coupling" begin
        E1 = 10.0
        ν1 = 0.33
        r_outer1 = 0.5
        r_fiber = 1.0
        data1 = LayerData(E1, ν1, r_outer1)
        var_indices = 10:17

        layer1 = PlaneLayer(data1, r_fiber, var_indices)

        E2 = 8.0
        ν2 = 0.4
        r_outer2 = 0.7
        data2 = LayerData(E2, ν2, r_outer2)

        ref_form1 = similar(layer1.ϕ)
        
        G1 = E1/(2*(1+ν1))
        κ1 = 3 - 4*ν1
        ϕ1 = layer1.ϕ
        z_bar_Φ1 = layer1.outer_z_bar_Φ
        bar_ψ1 = layer1.outer_bar_ψ
        add!(ref_form1, ϕ1,       κ1/G1+0.0im)
        add!(ref_form1, z_bar_Φ1, -1.0/G1+0.0im)
        add!(ref_form1, bar_ψ1,   -1.0/G1+0.0im)
        
        layer2 = PlaneLayer(data2, layer1)

        ref_form2 = similar(layer2.ϕ)

        G2 = E2/(2*(1+ν2))
        κ2 = 3 - 4*ν2
        ϕ2 = layer2.ϕ
        z_bar_Φ2 = layer2.inner_z_bar_Φ
        bar_ψ2 = layer2.inner_bar_ψ
        add!(ref_form2, ϕ2,       κ2/G2+0.0im)
        add!(ref_form2, z_bar_Φ2, -1.0/G2+0.0im)
        add!(ref_form2, bar_ψ2,   -1.0/G2+0.0im)

        @test my_isapprox(ref_form1, ref_form2, atol=1e-14)

        ref_form2_outer = similar(layer2.ϕ)

        z_bar_Φ2 = layer2.outer_z_bar_Φ
        bar_ψ2 = layer2.outer_bar_ψ
        add!(ref_form2_outer, ϕ2,       κ2/G2+0.0im)
        add!(ref_form2_outer, z_bar_Φ2, -1.0/G2+0.0im)
        add!(ref_form2_outer, bar_ψ2,   -1.0/G2+0.0im)

        E3 = 12.0
        ν3 = 0.1
        r_outer3 = 1.0
        data3 = LayerData(E3, ν3, r_outer3)

        layer3 = PlaneLayer(data3, layer2)

        ref_form3 = similar(layer3.ϕ)

        G3 = E3/(2*(1+ν3))
        κ3 = 3 - 4*ν3
        ϕ3 = layer3.ϕ
        z_bar_Φ3 = layer3.inner_z_bar_Φ
        bar_ψ3 = layer3.inner_bar_ψ
        add!(ref_form3, ϕ3,       κ3/G3+0.0im)
        add!(ref_form3, z_bar_Φ3, -1.0/G3+0.0im)
        add!(ref_form3, bar_ψ3,   -1.0/G3+0.0im)

        @test my_isapprox(ref_form2_outer, ref_form3, atol=1e-14)
    end
end

function test2()
    @testset "Fiber construction" begin
        layers_data = Vector(
            [
                LayerData(10.0, 0.33, 0.5),
                LayerData(8.0, 0.4, 0.7),
                LayerData(12.0, 0.1, 1.0)
            ]
        )
        var_indices = 10:21
        fiber = PlaneFiber(layers_data, var_indices)
        
        #########################################
        # Проверяем совпадение индексов
        flag = true
        flag = flag && coincede_indices(fiber[1].ϕ, fiber[1].inner_z_bar_Φ)
        flag = flag && coincede_indices(fiber[1].ϕ, fiber[1].outer_z_bar_Φ)
        flag = flag && coincede_indices(fiber[1].ϕ, fiber[1].inner_bar_ψ)
        flag = flag && coincede_indices(fiber[1].ϕ, fiber[1].outer_bar_ψ)
        for l in firstindex(fiber)+1 : lastindex(fiber)
            flag = flag && coincede_indices(fiber[l-1].ϕ,             fiber[l].ϕ)
            flag = flag && coincede_indices(fiber[l-1].inner_z_bar_Φ, fiber[l].inner_z_bar_Φ)
            flag = flag && coincede_indices(fiber[l-1].inner_bar_ψ,   fiber[l].inner_bar_ψ)
            flag = flag && coincede_indices(fiber[l-1].outer_z_bar_Φ, fiber[l].outer_z_bar_Φ)
            flag = flag && coincede_indices(fiber[l-1].outer_bar_ψ,   fiber[l].outer_bar_ψ)
        end      
        @test flag

        #########################################

        for l in firstindex(fiber)+1 : lastindex(fiber)
            layer1 = fiber[l-1]
            E1 = layer1.E
            ν1 = layer1.ν
            G1 = E1/(2*(1+ν1))
            κ1 = 3 - 4*ν1
            ϕ1 = layer1.ϕ
            z_bar_Φ1 = layer1.outer_z_bar_Φ
            bar_ψ1 = layer1.outer_bar_ψ
            
            displ_form1 = similar(ϕ1)
            add!(displ_form1, ϕ1,       κ1/G1+0.0im)
            add!(displ_form1, z_bar_Φ1, -1.0/G1+0.0im)
            add!(displ_form1, bar_ψ1,   -1.0/G1+0.0im)

            force_form1 = similar(ϕ1)
            add!(force_form1, ϕ1,       1.0+0.0im)
            add!(force_form1, z_bar_Φ1, 1.0+0.0im)
            add!(force_form1, bar_ψ1,   1.0+0.0im)


            layer2 = fiber[l]
            E2 = layer2.E
            ν2 = layer2.ν
            G2 = E2/(2*(1+ν2))
            κ2 = 3 - 4*ν2
            ϕ2 = layer2.ϕ
            z_bar_Φ2 = layer2.inner_z_bar_Φ
            bar_ψ2 = layer2.inner_bar_ψ
            
            displ_form2 = similar(ϕ2)
            add!(displ_form2, ϕ2,       κ2/G2+0.0im)
            add!(displ_form2, z_bar_Φ2, -1.0/G2+0.0im)
            add!(displ_form2, bar_ψ2,   -1.0/G2+0.0im)

            force_form2 = similar(ϕ2)
            add!(force_form2, ϕ2,       1.0+0.0im)
            add!(force_form2, z_bar_Φ2, 1.0+0.0im)
            add!(force_form2, bar_ψ2,   1.0+0.0im)

            @test my_isapprox(displ_form1, displ_form2, atol=1e-14)
            @test my_isapprox(force_form1, force_form2, atol=1e-14)
        end
    end
end

function test3() 
    @testset "Fiber contour forces and displacements" begin
        layers_data = Vector(
            [
                LayerData(10.0, 0.33, 0.5),
                LayerData(8.0, 0.4, 0.7),
                LayerData(12.0, 0.1, 1.0)
            ]
        )
        var_indices = 10:21
        fiber = PlaneFiber(layers_data, var_indices)

        layer = fiber[end]
        E = layer.E
        ν = layer.ν
        G = E/(2*(1+ν))
        κ = 3 - 4*ν
        ϕ = layer.ϕ
        z_bar_Φ = layer.outer_z_bar_Φ
        bar_ψ = layer.outer_bar_ψ
        
        displ_form = similar(ϕ)
        add!(displ_form, ϕ,       κ/(2G))
        add!(displ_form, z_bar_Φ, -1.0/(2G))
        add!(displ_form, bar_ψ,   -1.0/(2G))
        output = similar(ϕ[end])

        flag = true
        for v in var_indices
            displacements_series!(output, fiber, v)
            flag = flag && my_isapprox(displ_form[v], output, atol=1e-10)
        end
        @test flag

        force_form = similar(ϕ)
        add!(force_form, ϕ,       1.0)
        add!(force_form, z_bar_Φ, 1.0)
        add!(force_form, bar_ψ,   1.0)
        
        flag = true
        for v in var_indices
            forces_series!(output, fiber, v)
            flag = flag && my_isapprox(force_form[v], output, atol=1e-10)
        end
        @test flag
    end 
end

function test4()
    @testset "Cohesive" begin
        lattice = Lattice(1.0+0.0im, exp(1.0im))
        praecursor = EllipticPraecursor(lattice, 20)
        inclusions = [
            InclusionData(1//5+1//5im, 0.1, 10),
            InclusionData(4//5+4//5im, 0.1, 10),
            InclusionData(1//2+1//2im, 0.1, 10)
        ]
        E = 10.0
        ν = 0.3
        first_index = 43 
        cohesive = PlaneCohesive(E, ν, inclusions, first_index, praecursor)
    end

end

function test()
    test1()
    test2()
    test3()
    test4()
end

end # module TestPlaneProblems
using .TestPlaneProblems
TestPlaneProblems.test()