module TestPlaneProblems

using Test
using DelimitedFiles
using OffsetArrays

using VanFoFy.Testing: my_isapprox, coincede_indices
using VanFoFy.PlaneProblems: PlaneLayer, add!, PlaneFiber
using VanFoFy.PlaneProblems: forces_coupling!, displacements_coupling!
using VanFoFy.Input: LayerData, InclusionData, CellData, CohesiveData, FiberData
using VanFoFy.FunctionalTerms: EllipticPraecursor
using VanFoFy.Types: Lattice, BoundedVector
using VanFoFy.PlaneProblems: PlaneCohesive, PlaneProblem, eachpower
using VanFoFy.FunctionalTerms: ZTerm, ConstTerm, WeierstrassTerm, QSpecialTerm
using VanFoFy.FunctionalTerms: conjugate

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
        n_layers = size(fiber.layers, 1)
        flag = true
        flag = flag && coincede_indices(fiber.layers[1].ϕ, fiber.layers[1].inner_z_bar_Φ)
        flag = flag && coincede_indices(fiber.layers[1].ϕ, fiber.layers[1].outer_z_bar_Φ)
        flag = flag && coincede_indices(fiber.layers[1].ϕ, fiber.layers[1].inner_bar_ψ)
        flag = flag && coincede_indices(fiber.layers[1].ϕ, fiber.layers[1].outer_bar_ψ)
        for l in 2 : n_layers
            flag = flag && coincede_indices(fiber.layers[l-1].ϕ,             fiber.layers[l].ϕ)
            flag = flag && coincede_indices(fiber.layers[l-1].inner_z_bar_Φ, fiber.layers[l].inner_z_bar_Φ)
            flag = flag && coincede_indices(fiber.layers[l-1].inner_bar_ψ,   fiber.layers[l].inner_bar_ψ)
            flag = flag && coincede_indices(fiber.layers[l-1].outer_z_bar_Φ, fiber.layers[l].outer_z_bar_Φ)
            flag = flag && coincede_indices(fiber.layers[l-1].outer_bar_ψ,   fiber.layers[l].outer_bar_ψ)
        end      
        @test flag

        #########################################

        for l in 2 : n_layers
            layer1 = fiber.layers[l-1]
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


            layer2 = fiber.layers[l]
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

        layer = fiber.layers[end]
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
            displacements_coupling!(output, fiber, v)
            flag = flag && my_isapprox(displ_form[v], output, atol=1e-10)
        end
        @test flag

        force_form = similar(ϕ)
        add!(force_form, ϕ,       1.0)
        add!(force_form, z_bar_Φ, 1.0)
        add!(force_form, bar_ψ,   1.0)
        
        flag = true
        for v in var_indices
            forces_coupling!(output, fiber, v)
            flag = flag && my_isapprox(force_form[v], output, atol=1e-10)
        end
        @test flag
    end 

@testset "Multicoated fiber" begin
    ω1 = 1.0 + 0.0im
    ω3 = exp(1.0im)
    lattice = Lattice(ω1, ω3)
    praecursor = EllipticPraecursor(lattice, 6)

    cohesive_data = CohesiveData(1.0, 0.33)

    fibers_data = [
        FiberData(
            [13.5, 12.2, 11.8, 10.0],
            [0.16, 0.24, 0.26, 0.2],
            [0.13, 0.28, 0.34, 0.5],
            0.0,
            complex(1//3,1//3),
            6
        ),
        FiberData(
            [8.6, 11.0],
            [0.37, 0.3],
            [0.16, 0.2],
            0.0,
            complex(1//3, 5//7),
            5
        ),
        FiberData(
            [12.3],
            [0.4],
            [0.3],
            0.0,
            complex(5//7, 2//3),
            4
        )
    ]

    cell_data = CellData(
        2.0,
        2.0,
        1.0,
        fibers_data,
        cohesive_data
    )

    plane_problem = PlaneProblem(cell_data, praecursor)

    fiber = plane_problem.fibers[1]

    ref_vector = [
        1.0, 
        0.98524492234169658289, 
        0.97918549041852409598, 
        0.83967230348918997507]
    ind = firstindex(fiber.layers[1])
    comp_vector = [fiber.layers[k].ϕ[ind][0] for k in 1 : 4]

    # @test my_isapprox(ref_vector, comp_vector)

    ref_vector = [
        1.0, 
        1.0361001068980695283, 
        1.0464517928694767601, 
        0.93295839017634774049]
    ind = firstindex(fiber.layers[1]) + 2
    comp_vector = [fiber.layers[k].ϕ[ind][1] for k in 1 : 4]
    # @test my_isapprox(ref_vector, comp_vector)

    ref_vector = [
        1.0000000000000000000, 
        0.98524492234169658289, 
        0.97918549041852409598, 
        0.83967230348918997507]
    ind = firstindex(fiber.layers[1]) + 4
    comp_vector = [fiber.layers[k].ϕ[ind][2] for k in 1 : 4]
    # @test my_isapprox(ref_vector, comp_vector)

    function read_test_data(name::String)
        f_re = readdlm("test_data/" * name * "_re.txt", ',', Float64)
        f_im = readdlm("test_data/" * name * "_im.txt", ',', Float64)
        OffsetArray(
            f_re .+ f_im .* 1.0im,
            eachpower(fiber),
            eachindex(fiber)
        )
    end

    @test my_isapprox(fiber.layers[1].ϕ, read_test_data("phi_1_1"))
    @test my_isapprox(fiber.layers[1].outer_z_bar_Φ, read_test_data("z_bar_Phi_1_1"))
    @test my_isapprox(fiber.layers[1].outer_bar_ψ, read_test_data("bar_psi_1_1"))

    @test my_isapprox(fiber.layers[2].ϕ, read_test_data("phi_1_2"))
    @test my_isapprox(fiber.layers[2].outer_z_bar_Φ, read_test_data("z_bar_Phi_1_2"))
    @test my_isapprox(fiber.layers[2].outer_bar_ψ, read_test_data("bar_psi_1_2"))

    r1 = fiber.layers[1].r_outer
    r2 = fiber.layers[2].r_outer

    ϕ1       = fiber.layers[1].ϕ
    z_bar_Φ1 = fiber.layers[1].outer_z_bar_Φ
    bar_ψ1   = fiber.layers[1].outer_bar_ψ
    ϕ2       = fiber.layers[2].ϕ
    z_bar_Φ2 = fiber.layers[2].inner_z_bar_Φ

    bar_ψ2 = similar(fiber.layers[1].outer_bar_ψ)
    add!(bar_ψ2, ϕ1,       +1.0)
    add!(bar_ψ2, z_bar_Φ1, +1.0)
    add!(bar_ψ2, bar_ψ1,   +1.0)
    add!(bar_ψ2, ϕ2,       -1.0)
    add!(bar_ψ2, z_bar_Φ2, -1.0)
    bar_ψ2 = conjugate(bar_ψ2, r1)
    bar_ψ2 = conjugate(bar_ψ2, r2)

    @test my_isapprox(bar_ψ2, read_test_data("psi_test"))
    @test my_isapprox(bar_ψ2, read_test_data("bar_psi_1_2"))
    @test my_isapprox(fiber.layers[2].outer_bar_ψ, bar_ψ2)
end
end

function test4()
@testset "Cohesive" begin
    ω1 = 1.0 + 0.0im
    ω3 = exp(1.0im)
    lattice = Lattice(ω1, ω3)
    praecursor = EllipticPraecursor(lattice, 6)

    cohesive_data = CohesiveData(1.0, 0.33)

    fibers_data = [
        FiberData(
            [13.5, 12.2, 11.8, 10.0],
            [0.16, 0.24, 0.26, 0.2],
            [0.13, 0.28, 0.34, 0.5],
            0.0,
            complex(1//3,1//3),
            6
        ),
        FiberData(
            [8.6, 11.0],
            [0.37, 0.3],
            [0.16, 0.2],
            0.0,
            complex(1//3, 5//7),
            5
        ),
        FiberData(
            [12.3],
            [0.4],
            [0.3],
            0.0,
            complex(5//7, 2//3),
            4
        )
    ]

    cell_data = CellData(
        2.0,
        2.0,
        1.0,
        fibers_data,
        cohesive_data
    )

    plane_problem = PlaneProblem(cell_data, praecursor)

    cohesive = plane_problem.cohesive

    buffer = BoundedVector{ComplexF64}(-6 : 8)

    displacements_coupling!(
        buffer,
        cohesive,
        complex(1//3,1//3),
        0.5,
        3
    )

    ref_vector = OffsetVector(
        [
            0.000013137314474721333444 + 0.000060922089524345094324im, 
            0.000055595971328228494239 + 0.000064370198816722466457im, 
            -0.022717409856195618784 - 0.008509615513728190084im, 
            0.112461165037357313512 - 0.016030950306181172837im, 
            0.019326246046475602930 - 0.014437136718637749722im, 
            -2.2301102991334245118 - 0.0093731673949354055im, 
            0.0im, 
            0.0im, 
            -2.7314694224061679684 - 1.4922085463926697546im, 
            -1.33240223248528244504 - 0.00524897374116382702im, 
            0.0im, 
            0.0034481454838687932395 + 0.0004915212187426809499im, 
            0.0im, 
            -4.1777283301521625855e-6 + 4.8370627725249553200e-6im, 
            0.0im
        ],
        -6:8
    )

    @test my_isapprox(ref_vector, buffer, atol=1e-14)
end
end

function test5()
    ω1 = 1.0 + 0.0im
    ω3 = exp(1.0im)
    lattice = Lattice(ω1, ω3)
    praecursor = EllipticPraecursor(lattice, 6)

    cohesive_data = CohesiveData(1.0, 0.33)

    fibers_data = [
        FiberData(
            [13.5, 12.2, 11.8, 10.0],
            [0.16, 0.24, 0.26, 0.2],
            [0.13, 0.28, 0.34, 0.5],
            0.0,
            complex(1//3,1//3),
            6
        ),
        FiberData(
            [8.6, 11.0],
            [0.37, 0.3],
            [0.16, 0.2],
            0.0,
            complex(1//3, 5//7),
            5
        ),
        FiberData(
            [12.3],
            [0.4],
            [0.3],
            0.0,
            complex(5//7, 2//3),
            4
        )
    ]

    cell_data = CellData(
        2.0,
        2.0,
        1.0,
        fibers_data,
        cohesive_data
    )

    plane_problem = PlaneProblem(cell_data, praecursor)

    ref_matrix = readdlm("test_data/matrix.txt", ',', Float64)

    # @test my_isapprox(ref_matrix[:,1], plane_problem.matrix[1:end-4, 1], atol=1e-14)
    # @test my_isapprox(ref_matrix[:,2], plane_problem.matrix[1:end-4, 2], atol=1e-14)
    # @test my_isapprox(ref_matrix[:,3], plane_problem.matrix[1:end-4, 3], atol=1e-14)

    flag = true
    errors = 0
    for i in axes(ref_matrix, 2)
        ff = my_isapprox(ref_matrix[1:end-4,i], plane_problem.matrix[1:end-4, i], atol=1e-12)
        if !ff
            println("Column: ", i)
            errors += 1
        end
        flag = flag && ff
    end
    @test flag
end

function test()
    test1()
    test2()
    test3()
    test4()
    test5()
end

end # module TestPlaneProblems
using .TestPlaneProblems
TestPlaneProblems.test()