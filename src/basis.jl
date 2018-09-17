"""
Simple flat basis (all symmetry blocks) intended for iteration.
"""
struct FlatBasis{N}
    "Local basis truncation."
    l_max::Int
    "Many-body basis truncation."
    l_total_max::Int
end

mutable struct FlatBasisIterState
    v::Vector{Int}
    l_total::Int
end

function FlatBasisIterState(::FlatBasis{N}) where {N}
    v = zeros(Int, 2N)
    # m = -1
    v[2] -= 1

    l_total = 0

    FlatBasisIterState(v, l_total)
end

function Base.iterate(fb::FlatBasis{N}, state::FlatBasisIterState=FlatBasisIterState(fb)) where {N}
    for i in 1:N
        # m += 1
        state.v[2i] += 1

        # m > l
        if state.v[2i] > state.v[2i-1]
            state.l_total += 1
            # l += 1
            state.v[2i-1] += 1

            # l > l_max || l_total > l_total_max
            if state.v[2i-1] > fb.l_max || state.l_total > fb.l_total_max
                state.l_total -= state.v[2i-1]
                # l = 0
                state.v[2i-1] = 0
                # m = 0
                state.v[2i] = 0

                i == N && return nothing

                continue
            else
                # m = -l
                state.v[2i] = -state.v[2i-1]
            end
        end

        break
    end

    state.v, state
end

Base.eltype(::Type{FlatBasis}) = Vector{Int}


"""
Symmetry block label.
"""
struct BlockLabel
    "Total l parity."
    lp::Int
    "Total m value."
    m::Int
end


"""
Symmetry block of a basis for `N` rotors.
"""
struct Block{N}
    "Basis vectors, each of the form [l_1, m_1, ..., l_N, m_N]."
    vectors::Vector{Vector{Int}}
    "Mapping from basis vectors to their indices."
    lookup::Dict{Vector{Int},Int}

    "Number of basis vectors."
    size::Int
end

"""
    Block(N::Int, vectors::Vector{Vector{Int}}, lookup::Dict{Vector{Int},Int})

Create a block of `vectors` with corresponding `lookup` table for `N` rotors.
"""
function Block(N::Int, vectors::Vector{Vector{Int}}, lookup::Dict{Vector{Int},Int})
    # For non-empty blocks, make sure N is correct.
    length(vectors) == 0 || length(vectors[1]) == 2N || throw(DomainError(N, "Invalid number of rotors."))

    Block{N}(vectors, lookup, length(vectors))
end


abstract type AbstractBasis{N} end


"""
Truncated basis for `N` rotors, containing a single symmetry block.
"""
struct SingleBlockBasis{N} <: AbstractBasis{N}
    "Local basis truncation."
    l_max::Int
    "Many-body basis truncation."
    l_total_max::Int

    "Symmetry label for the block."
    label::BlockLabel
    "The block of basis vectors."
    block::Block{N}

    "Total number of basis vectors."
    size::Int
end

"""
    SingleBlockBasis(N::Int, l_max::Int, l_total_max::Int, sym_lp::Int, sym_m::Int)

Generate a truncated `N`-body basis of spherical harmonics.

The local bases are truncated at `l_max`. The many-body basis is truncated so
that the sum of `l` values does not exceed `l_total_max`.

Only states with total l parity `sym_lp` and total m value `sym_m` are
included.
"""
function SingleBlockBasis(N::Int, l_max::Int, l_total_max::Int, sym_lp::Int, sym_m::Int)
    N >= 1 || throw(DomainError(N, "At least one rotor."))
    l_max >= 0 || throw(DomainError(l_max, "Non-positive l_max."))
    l_total_max >= 0 || throw(DomainError(l_total_max, "Non-positive l_total_max."))

    vectors = Vector{Int}[]
    lookup = Dict{Vector{Int},Int}()

    for v in FlatBasis{N}(l_max, l_total_max)
        lp_total = sum(v[2i-1] for i in 1:N) % 2
        m_total = sum(v[2i] for i in 1:N)

        lp_total == sym_lp || continue
        m_total == sym_m || continue

        vv = copy(v)
        push!(vectors, vv)
        lookup[vv] = length(vectors)
    end

    label = BlockLabel(sym_lp, sym_m)
    block = Block(N, vectors, lookup)

    SingleBlockBasis{N}(l_max, l_total_max, label, block, block.size)
end


"""
Truncated basis for `N` rotors, containing multiple symmetry blocks.
"""
struct MultiBlockBasis{N} <: AbstractBasis{N}
    "Local basis truncation."
    l_max::Int
    "Many-body basis truncation."
    l_total_max::Int

    "Blocks of basis vectors."
    blocks::Dict{BlockLabel,Block{N}}

    "Total number of basis vectors."
    size::Int
end

"""
    MultiBlockBasis(N::Int, l_max::Int, l_total_max::Int)

Generate a truncated `N`-body basis of spherical harmonics.

The local bases are truncated at `l_max`. The many-body basis is truncated so
that the sum of `l` values does not exceed `l_total_max`.
"""
function MultiBlockBasis(N::Int, l_max::Int, l_total_max::Int)
    N >= 1 || throw(DomainError(N, "At least one rotor."))
    l_max >= 0 || throw(DomainError(l_max, "Non-positive l_max."))
    l_total_max >= 0 || throw(DomainError(l_total_max, "Non-positive l_total_max."))

    # For each site, the magnitude of m is bounded by l, so (by the triangle
    # inequality) the magnitude of the sum of the m values can't exceed the sum
    # of the l values, which is in turn bounded by l_total_max.
    vectors = Dict(((lp, m), Vector{Int}[]) for lp in 0:1 for m in -l_total_max:l_total_max)
    lookups = Dict(((lp, m), Dict{Vector{Int},Int}()) for lp in 0:1 for m in -l_total_max:l_total_max)

    for v in FlatBasis{N}(l_max, l_total_max)
        lp_total = sum(v[2i-1] for i in 1:N) % 2
        m_total = sum(v[2i] for i in 1:N)

        vv = copy(v)
        push!(vectors[(lp_total, m_total)], vv)
        lookups[(lp_total, m_total)][vv] = length(vectors[(lp_total, m_total)])
    end

    blocks = Dict{BlockLabel,Block{N}}()
    size = 0
    for (lp, m) in keys(vectors)
        label = BlockLabel(lp, m)
        blocks[label] = Block(N, vectors[(lp, m)], lookups[(lp, m)])
        size += blocks[label].size
    end

    MultiBlockBasis{N}(l_max, l_total_max, blocks, size)
end
