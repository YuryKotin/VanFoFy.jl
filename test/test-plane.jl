module TestPlaneProblems

using Test

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

        for v in eachindex(ref_form1)
            for p in eachindex(ref_form1[v])
                @test ref_form1[v][p] ≈ ref_form2[v][p]
            end
        end

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

        for v in eachindex(ref_form3)
            for p in eachindex(ref_form3[v])
                if ref_form2_outer[v][p] == NaN + NaN*im
                    println("2", v, p)
                end
                if ref_form3[v][p] == NaN + NaN*im
                    println("3", v, p)
                end
                @test ref_form2_outer[v][p] ≈ ref_form3[v][p] atol=1e-14
            end
        end
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
        
        for v in var_indices
            @test firstindex(fiber[1].ϕ[v]) == firstindex(fiber[1].inner_z_bar_Φ[v])   
            @test firstindex(fiber[1].ϕ[v]) == firstindex(fiber[1].outer_z_bar_Φ[v])
            @test firstindex(fiber[1].ϕ[v]) == firstindex(fiber[1].inner_bar_ψ[v])
            @test firstindex(fiber[1].ϕ[v]) == firstindex(fiber[1].outer_bar_ψ[v])      
        end
        for l in firstindex(fiber)+1 : lastindex(fiber)
            layer1 = fiber[l-1]
            layer2 = fiber[l]
            for v in var_indices
                @test firstindex(layer1.ϕ[v])             == firstindex(layer2.ϕ[v])
                @test firstindex(layer1.inner_z_bar_Φ[v]) == firstindex(layer2.inner_z_bar_Φ[v])
                @test firstindex(layer1.inner_bar_ψ[v])   == firstindex(layer2.inner_bar_ψ[v])
                @test firstindex(layer1.outer_z_bar_Φ[v]) == firstindex(layer2.outer_z_bar_Φ[v])
                @test firstindex(layer1.outer_bar_ψ[v])   == firstindex(layer2.outer_bar_ψ[v])
                
                @test lastindex(layer1.ϕ[v])             == lastindex(layer2.ϕ[v])
                @test lastindex(layer1.inner_z_bar_Φ[v]) == lastindex(layer2.inner_z_bar_Φ[v])
                @test lastindex(layer1.inner_bar_ψ[v])   == lastindex(layer2.inner_bar_ψ[v])
                @test lastindex(layer1.outer_z_bar_Φ[v]) == lastindex(layer2.outer_z_bar_Φ[v])
                @test lastindex(layer1.outer_bar_ψ[v])   == lastindex(layer2.outer_bar_ψ[v])
            end
        end      

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

            for v in eachindex(displ_form1)
                for p in eachindex(displ_form1[v])
                    @test displ_form1[v][p] ≈ displ_form2[v][p] atol=1e-14
                    @test force_form1[v][p] ≈ force_form2[v][p] atol=1e-14
                end
            end
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
        for v in var_indices
            displacements_series!(output, fiber, v)
            term = displ_form[v]
            for p in eachindex(output)
                @test term[p] ≈ output[p] atol=1e-10
            end
        end

        force_form = similar(ϕ)
        add!(force_form, ϕ,       1.0)
        add!(force_form, z_bar_Φ, 1.0)
        add!(force_form, bar_ψ,   1.0)
        for v in var_indices
            forces_series!(output, fiber, v)
            term = force_form[v]
            for p in eachindex(output)
                @test term[p] ≈ output[p] atol=1e-10
            end
        end
    end 
end

function test4()
    lattice = Lattice(1.0+0.0im, exp(1.0im))
    praecursor = EllipticPraecursor(lattice, 20)
    inclusions = [
        InclusionData(1//5+1//5im, 0.1, 10),
        InclusionData(4//5+4//5im, 0.1, 10),
        InclusionData(1//2+1//2im, 0.1, 10)
    ]
    E = 10.0
    ν = 0.3
    cohesive = PlaneCohesive(E, ν, inclusions, 43, praecursor)
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