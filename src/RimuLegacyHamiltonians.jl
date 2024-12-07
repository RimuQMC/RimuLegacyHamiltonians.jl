module RimuLegacyHamiltonians

using Rimu: Rimu, BitStringAddresses, SingleComponentFockAddress
using Rimu.Interfaces: AbstractHamiltonian, IsHermitian, AbstractOperator
using Rimu.BitStringAddresses: num_components, onr, near_uniform, time_reverse,
    print_address, AbstractFockAddress, BoseFS, OccupiedModeMap, CompositeFS
using Rimu.Hamiltonians: Hamiltonians, dimension, HubbardMom1D, starting_address,
    number_conserving_dimension, LOStructure, num_offdiagonals, diagonal_element,
    num_modes, num_occupied_modes, offdiagonals, get_offdiagonal, AbstractOffdiagonals,
    momentum_transfer_excitation, hopnextneighbour, HubbardReal1D, check_tr_address,
    num_singly_doubly_occupied_sites, DensityMatrixDiagonal, find_mode,
    number_conserving_dimension, Momentum, G2MomCorrelator

export BoseFS2C
export BoseHubbardMom1D2C, BoseHubbardReal1D2C, G2MomCorrelator

include("BoseFS2C.jl")
include("BoseHubbardMom1D2C.jl")
include("BoseHubbardReal1D2C.jl")
include("G2MomCorrelator.jl")

end
