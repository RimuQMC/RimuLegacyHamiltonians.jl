using RimuLegacyHamiltonians
using Rimu: num_modes, num_components, num_particles, CompositeFS, onr, BoseFS,
    near_uniform, starting_address, LOStructure, offdiagonals, IsHermitian, FermiFS, DVec,
    HubbardRealSpace, ExactDiagonalizationProblem
using Rimu.BitStringAddresses: parse_address
using Test

function exact_energy(ham)
    res = Rimu.solve(Rimu.ExactDiagonalizationProblem(ham))
    return res.values[1]
end


@testset "RimuLegacyHamiltonians.jl" begin
    @testset "BoseFS2C" begin
        fs1 = BoseFS2C((1,1,1,0,0,0,0), (1,1,1,0,5,0,0))
        @test num_modes(fs1) == 7
        @test num_components(fs1) == 2
        @test num_particles(fs1) == 11
        @test eval(Meta.parse(repr(fs1))) == fs1
        @test BoseFS2C(parse_address(sprint(show, fs1; context=:compact => true))) == fs1
        @test onr(fs1) == ([1, 1, 1, 0, 0, 0, 0], [1, 1, 1, 0, 5, 0, 0])

        fs2 = BoseFS2C(BoseFS((0,0,0,0,0,0,3)), BoseFS((0,2,1,0,5,0,0)))
        @test fs1 < fs2

        @test_throws MethodError BoseFS2C(BoseFS((1,1)), BoseFS((1,1,1)))
        @test BoseFS2C(CompositeFS(BoseFS((1,2)), BoseFS((3,1)))) == BoseFS2C((1,2), (3,1))
        @test CompositeFS(BoseFS2C(BoseFS((1,2)), BoseFS((3,1)))) ==
            CompositeFS(BoseFS((1,2)), BoseFS((3,1)))
    end

    @testset "2C model properties" begin
        flip(b) = BoseFS2C(b.bsb, b.bsa)
        addr1 = near_uniform(BoseFS2C{1,100,20})
        addr2 = near_uniform(BoseFS2C{100,1,20})

        for Hamiltonian in (BoseHubbardReal1D2C, BoseHubbardMom1D2C)
            @testset "$Hamiltonian" begin
                H1 = BoseHubbardReal1D2C(addr1; ta=1.0, tb=2.0, ua=0.5, ub=0.7, v=0.2)
                H2 = BoseHubbardReal1D2C(addr2; ta=2.0, tb=1.0, ua=0.7, ub=0.5, v=0.2)
                @test starting_address(H1) == addr1
                @test LOStructure(H1) == IsHermitian()

                hops1 = collect(offdiagonals(H1, addr1))
                hops2 = collect(offdiagonals(H2, addr2))
                sort!(hops1, by=a -> first(a).bsa)
                sort!(hops2, by=a -> first(a).bsb)

                addrs1 = first.(hops1)
                addrs2 = flip.(first.(hops2))
                values1 = last.(hops1)
                values2 = last.(hops1)
                @test addrs1 == addrs2
                @test values1 == values2

                @test eval(Meta.parse(repr(H1))) == H1
                @test eval(Meta.parse(repr(H2))) == H2
            end
        end
    end
    @testset "1D Bosons (2-component)" begin
        add1 = BoseFS2C(
            (1, 1, 1, 0, 0, 0),
            (1, 0, 0, 0, 0, 0),
        )
        H1 = BoseHubbardReal1D2C(add1, ua=2, v=3, tb=4)

        add2 = CompositeFS(
            BoseFS((1, 1, 1, 0, 0, 0)),
            BoseFS((1, 0, 0, 0, 0, 0)),
        )
        H2 = HubbardRealSpace(add2, t=[1, 4], u=[2 3; 3 0])

        add3 = CompositeFS(
            BoseFS((1, 1, 1, 0, 0, 0)),
            FermiFS((1, 0, 0, 0, 0, 0)),
        )
        H3 = HubbardRealSpace(add3, t=[1, 4], u=[2 3; 3 0])

        add4 = CompositeFS(
            BoseFS((1, 0, 0, 0, 0, 0)),
            BoseFS((1, 1, 1, 0, 0, 0)),
        )
        H4 = HubbardRealSpace(add4, t=[4, 1], u=[0 3; 3 2])

        add5 = CompositeFS(
            FermiFS((1, 0, 0, 0, 0, 0)),
            BoseFS((1, 1, 1, 0, 0, 0)),
        )
        H5 = HubbardRealSpace(add5, t=[4, 1], u=[0 3; 3 2])

        E1 = exact_energy(H1)
        E2 = exact_energy(H2)
        E3 = exact_energy(H3)
        E4 = exact_energy(H4)
        E5 = exact_energy(H5)

        @test E1 ≈ E2 rtol = 0.0001
        @test E2 ≈ E3 rtol = 0.0001
        @test E3 ≈ E4 rtol = 0.0001
        @test E4 ≈ E5 rtol = 0.0001
    end

end
