module TestQSpecial

using Test, OffsetArrays
using VanFoFy.Ellipticals: QSpecial, raw_complex, Weierstrass

function test()
    ω1 = complex(1.0)
    ω3 = exp(1im)
    wei = Weierstrass(ω1, ω3, 11)
    Q = QSpecial(wei, 10)

    rz = complex(3//10, 7//10)
    z = raw_complex(wei, rz)

    @test Q[rz, 0]  ≈ -0.51663997136492256779 + 0.59933186734118426564im
    @test Q[rz, 1]  ≈ 0.03059293447909888428 + 1.19592812759260613831im
    @test Q[rz, 2]  ≈ 0.61792409039690798789 + 0.53104522549332222603im
    @test Q[rz, 3]  ≈ 0.21949109770876260028 - 0.23037200859738624081im
    @test Q[rz, 4]  ≈ -0.25667024508338398547 - 0.26175052836614570717im
    @test Q[rz, 5]  ≈ -0.34680442894981433621 + 0.03239450713556108613im
    @test Q[rz, 6]  ≈ -0.09690359568209910845 + 0.21090154612274475410im
    @test Q[rz, 7]  ≈ 0.12370320077810871562 + 0.08401792981864268650im
    @test Q[rz, 8]  ≈ 0.100108873864963227041 - 0.090530908012789471084im
    @test Q[rz, 9]  ≈ -0.027628709900863662741 - 0.107911340307848935272im
    @test Q[rz, 10] ≈ -0.076637122750858210907 - 0.012691979234158738341im
end

end # module TestQSpecial
using .TestQSpecial
TestQSpecial.test()