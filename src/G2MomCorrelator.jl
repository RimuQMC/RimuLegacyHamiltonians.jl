
import Rimu.Hamiltonians: num_offdiagonals, diagonal_element, get_offdiagonal,
    G2MomCorrelator

"""
    G2MomCorrelator(d::Int, c=:cross) <: AbstractOperator{ComplexF64}

Two-body correlation operator representing the density-density correlation at distance `d`
of a two component system in a momentum-space Fock-state basis with addresses of type
[`BoseFS2C`](@ref). It returns a `Complex`
value.

Correlation across two components:
```math
\\hat{G}^{(2)}(d) = \\frac{1}{M}\\sum_{spqr=1}^M e^{-id(p-q)2π/M} a^†_{s} b^†_{p}  b_q a_r δ_{s+p,q+r}
```
Correlation within a single component:
```math
\\hat{G}^{(2)}(d) = \\frac{1}{M}\\sum_{spqr=1}^M e^{-id(p-q)2π/M} a^†_{s} a^†_{p}  a_q a_r δ_{s+p,q+r}
```

The diagonal element, where `(p-q)=0`, is
```math
\\frac{1}{M}\\sum_{k,p=1}^M a^†_{k} b^†_{p}  b_p a_k .
```

# Arguments
- `d::Integer`: the distance between two particles.
- `c`: possible instructions: `:cross`: default instruction, computing correlation between
  particles across two components;
  `:first`: computing correlation between particles within the first component;
  `:second`: computing correlation between particles within the second component.
  These are the only defined instructions, using anything else will produce errors.

# To use on a one-component system

For a system with only one component, e.g. with `BoseFS`, the second argument `c` is
irrelevant and can be any of the above instructions, one could simply skip this argument
and let it be the default value.

# See also

* [`BoseHubbardMom1D2C`](@ref)
* [`BoseFS2C`](@ref)
* [`Rimu.G2RealCorrelator`](@extref)
* [`Rimu.G2RealSpace`](@extref)
* [`Rimu.AbstractOperator`](@extref)
* [`Rimu.AllOverlaps`](@extref)
"""
function G2MomCorrelator(d, c)
    if c == :first
        return G2MomCorrelator{1}(d)
    elseif c == :second
        return G2MomCorrelator{2}(d)
    elseif c == :cross
        return G2MomCorrelator{3}(d)
    else
        throw(ArgumentError("Unknown instruction for G2MomCorrelator!"))
    end
end

function Rimu.Interfaces.allows_address_type(g2m::G2MomCorrelator, ::Type{A}) where {A<:BoseFS2C}
    return num_modes(A) > g2m.d
end

function Base.show(io::IO, g::G2MomCorrelator{C}) where {C}
    if C == 1
        print(io, "G2MomCorrelator($(g.d), :first)")
    elseif C == 2
        print(io, "G2MomCorrelator($(g.d), :second)")
    else # C == 3 is the default value and handled in Rimu.jl
        @error "G2MomCorrelator{$C} is not implemented!"
    end
end

num_offdiagonals(g::G2MomCorrelator{1},addr::BoseFS2C) = num_offdiagonals(g, addr.bsa)
num_offdiagonals(g::G2MomCorrelator{2},addr::BoseFS2C) = num_offdiagonals(g, addr.bsb)

function num_offdiagonals(g::G2MomCorrelator{3}, addr::BoseFS2C)
    m = num_modes(addr)
    sa = num_occupied_modes(addr.bsa)
    sb = num_occupied_modes(addr.bsb)
    return sa*(m-1)*sb
end

diagonal_element(g::G2MomCorrelator{1}, addr::BoseFS2C) = diagonal_element(g, addr.bsa)
diagonal_element(g::G2MomCorrelator{2}, addr::BoseFS2C) = diagonal_element(g, addr.bsb)

function diagonal_element(g::G2MomCorrelator{3}, addr::BoseFS2C{NA,NB,M,AA,AB}) where {NA,NB,M,AA,AB}
    onrep_a = onr(addr.bsa)
    onrep_b = onr(addr.bsb)
    gd = 0
    for p in 1:M
        iszero(onrep_b[p]) && continue
        for k in 1:M
            gd += onrep_a[k] * onrep_b[p] # b†_p b_p a†_k a_k
        end
    end
    return ComplexF64(gd / M)
end

function get_offdiagonal(g::G2MomCorrelator{1}, addr::A, chosen)::Tuple{A,ComplexF64} where {A<:BoseFS2C}
    new_bsa, elem = get_offdiagonal(g, addr.bsa, chosen)
    return A(new_bsa, addr.bsb), elem
end

function get_offdiagonal(g::G2MomCorrelator{2}, addr::A, chosen)::Tuple{A,ComplexF64} where {A<:BoseFS2C}
    new_bsb, elem = get_offdiagonal(g, addr.bsb, chosen)
    return A(addr.bsa, new_bsb), elem
end

function get_offdiagonal(
    g::G2MomCorrelator{3},
    addr::A,
    chosen,
    ma=OccupiedModeMap(addr.bsa),
    mb=OccupiedModeMap(addr.bsb),
)::Tuple{A,ComplexF64} where {A<:BoseFS2C}

    m = num_modes(addr)
    new_bsa, new_bsb, gamma, _, _, Δp = momentum_transfer_excitation(
        addr.bsa, addr.bsb, chosen, ma, mb
    )
    gd = exp(-im * g.d * Δp * 2π / m) * gamma
    return A(new_bsa, new_bsb), ComplexF64(gd / m)
end
