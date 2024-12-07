using RimuLegacyHamiltonians
using Rimu
using Rimu.BitStringAddresses: parse_address
using Rimu.InterfaceTests # requires Rimu v0.14.0
using Test

function exact_energy(ham)
    res = Rimu.solve(Rimu.ExactDiagonalizationProblem(ham))
    return res.values[1]
end

using RimuLegacyHamiltonians: bose_hubbard_2c_interaction
@testset "RimuLegacyHamiltonians.jl" begin
    @testset "BoseFS2C 1" begin
        bfs2c = BoseFS2C(BoseFS((1, 2, 0, 4)), BoseFS((4, 0, 3, 1)))
        @test typeof(bfs2c) <: BoseFS2C{7,8,4}
        @test num_occupied_modes(bfs2c.bsa) == 3
        @test num_occupied_modes(bfs2c.bsb) == 3
        @test onr(bfs2c.bsa) == [1, 2, 0, 4]
        @test onr(bfs2c.bsb) == [4, 0, 3, 1]
        @test bose_hubbard_2c_interaction(bfs2c) == 8 # n_a*n_b over all sites
    end

    @testset "BoseFS2C 2" begin
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

        for address in (BoseFS2C((0, 1, 2, 3, 0), (1, 2, 3, 4, 5)), FermiFS2C((1, 0, 0, 1), (0, 0, 1, 0)))
            @test diagonal_element(Momentum(1), address) + diagonal_element(Momentum(2), address) ≡
                  diagonal_element(Momentum(0), address)
        end
        for address in (
            CompositeFS(BoseFS((1, 2, 3, 4, 5)), BoseFS((5, 4, 3, 2, 1))),
            BoseFS2C((1, 2, 3, 4, 5), (5, 4, 3, 2, 1))
        )
            for i in 1:5
                @test diagonal_element(DensityMatrixDiagonal(i, component=1), address) == i
                @test diagonal_element(DensityMatrixDiagonal(i, component=2), address) == 6 - i
                @test diagonal_element(DensityMatrixDiagonal(i), address) == 6
            end
        end
        @test_throws ArgumentError TimeReversalSymmetry(BoseHubbardMom1D2C(BoseFS2C((1, 1), (2, 1))))
        @testset "2-particle BoseHubbardMom1D2C" begin
            ham = BoseHubbardMom1D2C(BoseFS2C((0, 1, 1), (1, 0, 1)))
            even = TimeReversalSymmetry(ham)
            odd = TimeReversalSymmetry(ham; even=false)

            h_eigs = eigvals(Matrix(ham))
            p_eigs = sort!(vcat(eigvals(Matrix(even)), eigvals(Matrix(odd))))

            @test starting_address(even) == time_reverse(starting_address(ham))
            @test h_eigs ≈ p_eigs

            @test issymmetric(Matrix(odd))
            @test issymmetric(Matrix(even))
            @test LOStructure(odd) isa IsHermitian
        end

        for address in (
            BoseFS2C((1, 2, 3), (0, 1, 0)),
            CompositeFS(BoseFS((1, 2, 3)), FermiFS((0, 1, 0)))
        )
            @test single_particle_density(address) == (1, 3, 3)
            @test single_particle_density(address; component=1) == (1, 2, 3)
            @test single_particle_density(address; component=2) == (0, 1, 0)
            @test single_particle_density(DVec(address => 1); component=2) == (0, 7, 0)
        end

    end
    @testset "TwoComponentBosonicHamiltonian" begin
        aIni2cReal = BoseFS2C(BoseFS((1, 1, 1, 1)), BoseFS((1, 1, 1, 1))) # real space two-component
        Ĥ2cReal = BoseHubbardReal1D2C(aIni2cReal; ua=6.0, ub=6.0, ta=1.0, tb=1.0, v=6.0)
        hamA = HubbardReal1D(BoseFS((1, 1, 1, 1)); u=6.0, t=1.0)
        hamB = HubbardReal1D(BoseFS((1, 1, 1, 1)); u=6.0)

        test_hamiltonian_interface(Ĥ2cReal)

        @test hamA == Ĥ2cReal.ha
        @test hamB == Ĥ2cReal.hb
        @test num_offdiagonals(Ĥ2cReal, aIni2cReal) == 16
        @test num_offdiagonals(Ĥ2cReal, aIni2cReal) == num_offdiagonals(Ĥ2cReal.ha, aIni2cReal.bsa) + num_offdiagonals(Ĥ2cReal.hb, aIni2cReal.bsb)
        @test dimension(Ĥ2cReal) == 1225
        @test dimension(Float64, Ĥ2cReal) == 1225.0

        hp2c = offdiagonals(Ĥ2cReal, aIni2cReal)
        @test length(hp2c) == 16
        @test hp2c[1][1] == BoseFS2C(BoseFS((0, 2, 1, 1)), BoseFS((1, 1, 1, 1)))
        @test hp2c[1][2] ≈ -1.4142135623730951
        @test diagonal_element(Ĥ2cReal, aIni2cReal) ≈ 24.0 # from the V term

        aIni2cMom = BoseFS2C(BoseFS((0, 4, 0, 0)), BoseFS((0, 4, 0, 0))) # momentum space two-component
        Ĥ2cMom = BoseHubbardMom1D2C(aIni2cMom; ua=6.0, ub=6.0, ta=1.0, tb=1.0, v=6.0)

        test_hamiltonian_interface(Ĥ2cMom)

        @test num_offdiagonals(Ĥ2cMom, aIni2cMom) == 9
        @test dimension(Ĥ2cMom) == 1225
        @test dimension(Float64, Ĥ2cMom) == 1225.0

        hp2cMom = offdiagonals(Ĥ2cMom, aIni2cMom)
        @test length(hp2cMom) == 9
        @test hp2cMom[1][1] == BoseFS2C(BoseFS((1, 2, 1, 0)), BoseFS((0, 4, 0, 0)))
        @test hp2cMom[1][2] ≈ 2.598076211353316

        smat2cReal, adds2cReal = ExactDiagonalization.build_sparse_matrix_from_LO(Ĥ2cReal, aIni2cReal)
        eig2cReal = eigen(Matrix(smat2cReal))
        smat2cMom, adds2cMom = ExactDiagonalization.build_sparse_matrix_from_LO(Ĥ2cMom, aIni2cMom)
        eig2cMom = eigen(Matrix(smat2cMom))
        @test eig2cReal.values[1] ≈ eig2cMom.values[1]
    end

    @testset "2C Hamiltonian interface and structure" begin
        fs2c = BoseFS2C((1, 0, 1), (0, 1, 0))
        hr = BoseHubbardReal1D2C(fs2c)
        test_hamiltonian_interface(hr)
        test_hamiltonian_structure(hr)

        hm = BoseHubbardMom1D2C(fs2c)
        test_hamiltonian_interface(hm)
        test_hamiltonian_structure(hm)
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
    @testset "G2MomCorrelator" begin
        # v0 is the exact ground state from BoseHubbardMom1D2C(aIni;ua=0,ub=0,v=0.1)
        bfs1 = BoseFS([0, 2, 0])
        bfs2 = BoseFS([0, 1, 0])
        aIni = BoseFS2C(bfs1, bfs2)
        v0 = DVec(
            BoseFS2C((0, 2, 0), (0, 1, 0)) => 0.9999389545691221,
            BoseFS2C((1, 1, 0), (0, 0, 1)) => -0.007812695959057453,
            BoseFS2C((0, 1, 1), (1, 0, 0)) => -0.007812695959057453,
            BoseFS2C((2, 0, 0), (1, 0, 0)) => 4.046694762039993e-5,
            BoseFS2C((0, 0, 2), (0, 0, 1)) => 4.046694762039993e-5,
            BoseFS2C((1, 0, 1), (0, 1, 0)) => 8.616127793651117e-5,
        )
        g0 = G2MomCorrelator(0)
        g1 = G2MomCorrelator(1)
        g2 = G2MomCorrelator(2)
        g3 = G2MomCorrelator(3)
        @test imag(dot(v0, g0, v0)) == 0 # should be strictly real
        @test abs(imag(dot(v0, g3, v0))) < 1e-10
        @test dot(v0, g0, v0) ≈ 0.65 rtol = 0.01
        @test dot(v0, g1, v0) ≈ 0.67 rtol = 0.01
        @test dot(v0, g2, v0) ≈ 0.67 rtol = 0.01
        @test dot(v0, g3, v0) ≈ 0.65 rtol = 0.01
        @test num_offdiagonals(g0, aIni) == 2

        # on first component
        g0f = G2MomCorrelator(0, :first)
        g1f = G2MomCorrelator(1, :first)
        @test imag(dot(v0, g0f, v0)) == 0 # should be strictly real
        @test dot(v0, g0f, v0) ≈ 1.33 rtol = 0.01
        @test dot(v0, g1f, v0) ≈ 1.33 + 7.08e-5im rtol = 0.01
        # on second component
        g0s = G2MomCorrelator(0, :second)
        g1s = G2MomCorrelator(1, :second)
        #@test_throws ErrorException("invalid ONR") get_offdiagonal(g0s,aIni,1) # should fail due to invalid ONR
        @test dot(v0, g0s, v0) ≈ 1 / 3
        @test dot(v0, g1s, v0) ≈ 1 / 3
        # test against BoseFS
        ham1 = HubbardMom1D(bfs1)
        ham2 = HubbardMom1D(bfs2)
        @test num_offdiagonals(g0f, aIni) == num_offdiagonals(ham1, bfs1)
        @test num_offdiagonals(g0s, aIni) == num_offdiagonals(ham2, bfs2)
        aIni = BoseFS2C(bfs2, bfs1) # flip bfs1 and bfs2
        @test get_offdiagonal(g0s, aIni, 1) == (BoseFS2C(BoseFS{1,3}((0, 1, 0)), BoseFS{2,3}((1, 0, 1))), 0.47140452079103173)
        # test on BoseFS
        @test diagonal_element(g0s, bfs1) == 4 / 3
        @test diagonal_element(g0s, bfs2) == 1 / 3
        test_operator_interface(G2MomCorrelator(3), BoseFS(1, 2, 0, 3, 0, 4, 0, 1))
        test_operator_interface(G2MomCorrelator(2), BoseFS2C(BoseFS{1,3}((0, 1, 0)), BoseFS{2,3}((1, 0, 1))))
    end

end
