"""
    BoseFS2C{NA,NB,M,AA,AB} <: AbstractFockAddress
    BoseFS2C(onr_a, onr_b)

Address type that constructed with two [`BoseFS{N,M,S}`](@ref). It represents a
Fock state with two components, e.g. two different species of bosons with particle
number `NA` from species S and particle number `NB` from species B. The number of
modes `M` is expected to be the same for both components.
"""
struct BoseFS2C{NA,NB,M,SA,SB,N} <: AbstractFockAddress{N,M}
    bsa::BoseFS{NA,M,SA}
    bsb::BoseFS{NB,M,SB}
end

function BoseFS2C(bsa::BoseFS{NA,M,SA}, bsb::BoseFS{NB,M,SB}) where {NA,NB,M,SA,SB}
    N = NA + NB
    return BoseFS2C{NA,NB,M,SA,SB,N}(bsa, bsb)
end
BoseFS2C(onr_a::Tuple, onr_b::Tuple) = BoseFS2C(BoseFS(onr_a),BoseFS(onr_b))

function Rimu.BitStringAddresses.print_address(io::IO, b::BoseFS2C; compact)
    if compact
        print_address(io, b.bsa; compact)
        print(io, " ⊗ ")
        print_address(io, b.bsb; compact)
    else
        print(io, "BoseFS2C(", b.bsa, ", ", b.bsb, ")")
    end
end

Rimu.BitStringAddresses.num_components(::Type{<:BoseFS2C}) = 2
Base.isless(a::T, b::T) where {T<:BoseFS2C} = isless((a.bsa, a.bsb), (b.bsa, b.bsb))
Rimu.BitStringAddresses.onr(b2::BoseFS2C) = (onr(b2.bsa), onr(b2.bsb))


function Rimu.BitStringAddresses.near_uniform(::Type{<:BoseFS2C{NA,NB,M}}) where {NA,NB,M}
    return BoseFS2C(near_uniform(BoseFS{NA,M}), near_uniform(BoseFS{NB,M}))
end

function Rimu.BitStringAddresses.time_reverse(c::BoseFS2C{NA,NA,M,S,S,N}) where {NA,M,S,N}
    return  BoseFS2C{NA,NA,M,S,S,N}(c.bsb, c.bsa)
end

function Rimu.Hamiltonians.check_tr_address(addr::BoseFS2C{NA,NB,M,SA,SB,N}) where {NA,NB,M,SA,SB,N}
    if NA ≠ NB || SA ≠ SB
        throw(ArgumentError("Two component address with equal particle numbers required for `TimeReversalSymmetry`."))
    end
end

function Rimu.Interfaces.diagonal_element(dmd::DensityMatrixDiagonal{0}, add::BoseFS2C)
    return float(find_mode(add.bsa, dmd.mode).occnum + find_mode(add.bsb, dmd.mode).occnum)
end
function Rimu.Interfaces.diagonal_element(dmd::DensityMatrixDiagonal{1}, add::BoseFS2C)
    comp = add.bsa
    return float(find_mode(comp, dmd.mode).occnum)
end
function Rimu.Interfaces.diagonal_element(dmd::DensityMatrixDiagonal{2}, add::BoseFS2C)
    comp = add.bsb
    return float(find_mode(comp, dmd.mode).occnum)
end

function Rimu.Hamiltonians.dimension(::Type{<:BoseFS2C{NA,NB,M}}) where {NA,NB,M}
    return dimension(BoseFS{NA,M}) * dimension(BoseFS{NB,M})
end
function Rimu.Hamiltonians.number_conserving_dimension(addr::BoseFS2C)
    return number_conserving_dimension(addr.bsa) * number_conserving_dimension(addr.bsb)
end

BoseFS2C(fs::CompositeFS{2}) = BoseFS2C(fs.components...)
Rimu.CompositeFS(fs::BoseFS2C) = CompositeFS(fs.bsa, fs.bsb)

function Rimu.Interfaces.diagonal_element(m::Momentum{0}, a::BoseFS2C)
    Rimu.Hamiltonians._momentum((a.bsa, a.bsb), m.fold)
end
function Rimu.Interfaces.diagonal_element(m::Momentum{1}, a::BoseFS2C)
    Rimu.Hamiltonians._momentum(a.bsa, m.fold)
end
function Rimu.Interfaces.diagonal_element(m::Momentum{2}, a::BoseFS2C)
    Rimu.Hamiltonians._momentum(a.bsb, m.fold)
end

function Rimu.single_particle_density(add::Union{BoseFS2C}; component=0)
    if component == 0
        return float.(Tuple(sum(onr(add))))
    else
        return float.(Tuple(onr(add)[component]))
    end
end
