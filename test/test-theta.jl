@testset "θ-functions" begin
    ω1 = complex(1.0)
    ω3 = exp(1im)
    q = nome(ω1=ω1, ω3=ω3)
    tol_check = 1e-15
    tol_compute = 1e-20


    z = 0.63256929128615824176 + 0.47122375149242207160im

    th = OffsetArray{ComplexF64}(undef, 1:4, 0:4)

    for n in 0:4
        for k in 1:4
            th[k, n] = theta(th_k=k, d_n=n, z=z, q=q, ϵ=tol_compute)
        end
    end

    #####

    th_ma = [ 0.46214180304899833678 + 0.65513824162092381551im,
               0.96467023166692638429 + 0.11945667130003175272im,
               1.13835201868809736642 + 0.08123099419582718546im,
               0.86149287144921489500 - 0.08153295492698743884im]

    @test th[:, 0] ≈ th_ma atol=tol_check

    #####

    d1_th_ma = [ 0.97727706225681589769 + 0.08049915546539055744im,
             -0.42045304182205492196 - 0.63824580596766239552im,
              0.14241824316243509356 - 0.38564309515712238239im,
             -0.14363258751912619862 + 0.38620176200689986572im]

    @test th[:, 1] ≈ d1_th_ma atol=tol_check

    #####

    d2_th_ma = [-0.54549882845031761232 - 0.68893198268462155795im,
             -0.93946558497506382639 - 0.19739257822094024267im,
             -0.55247742739635853599 - 0.32311219018650595886im,
              0.55495918519936412129 + 0.32794356188507020189im]

    @test th[:, 2] ≈ d2_th_ma atol=tol_check

    #####

    d3_th_ma = [-1.05294508933080125267 + 0.15318331373277377347im,
              0.17025898335751795341 + 0.53691780124847588798im,
             -0.5623867730174600491 + 1.5392204494203911042im,
              0.5818162827245192532 - 1.5481591190168449923im]

    @test th[:, 3] ≈ d3_th_ma atol=tol_check

    #####

    d4_th_ma = [1.2954660922120757234 + 0.9931820896285708616im,
             0.71273193831970297909 + 0.89906624330941167631im,
             2.1950195882908139481 + 1.2634597310003502748im,
            -2.2347277131390164682 - 1.3407616781773902616im]

    @test th[:, 4] ≈ d4_th_ma atol=tol_check

    #####

    nome_powers = precompute_nome_powers(q=q, ϵ=tol_compute, z0=ω3)
    th1, th2, dth1, dth2 = theta_1d2d(z=z, nome_powers=nome_powers)

    @test th1  ≈ theta(th_k=1, d_n=0, z=z, q=q, ϵ=tol_compute) atol=tol_check
    @test th2  ≈ theta(th_k=2, d_n=0, z=z, q=q, ϵ=tol_compute) atol=tol_check
    @test dth1 ≈ theta(th_k=1, d_n=1, z=z, q=q, ϵ=tol_compute) atol=tol_check
    @test dth2 ≈ theta(th_k=2, d_n=1, z=z, q=q, ϵ=tol_compute) atol=tol_check
end