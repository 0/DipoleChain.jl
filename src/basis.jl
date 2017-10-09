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
    v_prev::Vector{Int}
    v::Vector{Int}

    l_total::Int
    done::Bool
end

function Base.start{N}(::FlatBasis{N})
    v_prev = Array{Int}(2N)
    v = zeros(Int, 2N)

    l_total = 0
    done = false

    FlatBasisIterState(v_prev, v, l_total, done)
end

function Base.next{N}(fb::FlatBasis{N}, state::FlatBasisIterState)
    state.v_prev .= state.v

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

                i == N && (state.done = true)

                continue
            else
                # m = -l
                state.v[2i] = -state.v[2i-1]
            end
        end

        break
    end

    state.v_prev, state
end

Base.done(::FlatBasis, state::FlatBasisIterState) = state.done

Base.eltype(::Type{FlatBasis}) = Vector{Int}


"""
Truncated basis for `N` rotors.
"""
struct Basis{N}
    "Local basis truncation."
    l_max::Int
    "Many-body basis truncation."
    l_total_max::Int

    "Symmetry block total l parity."
    sym_lp::Int
    "Symmetry block total m value."
    sym_m::Int

    "Number of basis vectors."
    size::Int

    "Basis vectors, each of the form [l_1, m_1, ..., l_N, m_N]."
    vectors::Vector{Vector{Int}}
    "Mapping from basis vectors to their indices."
    lookup::Dict{Vector{Int},Int}
end

"""
    Basis(N::Int, l_max::Int, l_total_max::Int, sym_lp::Int, sym_m::Int)

Generate a truncated `N`-body basis of spherical harmonics.

The local bases are truncated at `l_max`. The many-body basis is truncated so
that the sum of `l` values does not exceed `l_total_max`.

Only states with total l parity `sym_lp` and total m value `sym_m` are
included.
"""
function Basis(N::Int, l_max::Int, l_total_max::Int, sym_lp::Int, sym_m::Int)
    # At least one rotor.
    N >= 1 || throw(DomainError())
    # Non-negative l.
    l_max >= 0 || throw(DomainError())
    l_total_max >= 0 || throw(DomainError())

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

    Basis{N}(l_max, l_total_max, sym_lp, sym_m, length(vectors), vectors, lookup)
end


"""
Truncated basis for a subsystem containing `NA` rotors.
"""
struct SubsystemBasis{NA}
    "Basis for the system of which this is a subsystem."
    parent::Basis

    "Symmetry block labels of the form (lp, m)."
    blocks::Vector{Tuple{Int,Int}}

    "Number of basis vectors in each block."
    sizes::Dict{Tuple{Int,Int},Int}
    "Blocks of basis vectors."
    vectors::Dict{Tuple{Int,Int},Vector{Vector{Int}}}
end

"""
    SubsystemBasis{N}(parent::Basis{N}, NA::Int)

Generate a basis for `NA` sites of the full `parent` basis.
"""
function SubsystemBasis{N}(parent::Basis{N}, NA::Int)
    # At least one rotor.
    NA >= 1 || throw(DomainError())
    # Not all the rotors.
    NA < N || throw(DomainError())

    l_max = parent.l_max
    l_total_max = parent.l_total_max

    # For each site, the magnitude of m is bounded by l, so (by the triangle
    # inequality) the magnitude of the sum of the m values can't exceed the sum
    # of the l values, which is in turn bounded by l_total_max.
    sizes = Dict(((lp, m), 0) for lp in 0:1 for m in -l_total_max:l_total_max)
    vectors = Dict(((lp, m), Vector{Int}[]) for lp in 0:1 for m in -l_total_max:l_total_max)

    for v in FlatBasis{NA}(l_max, l_total_max)
        lp_total = sum(v[2i-1] for i in 1:NA) % 2
        m_total = sum(v[2i] for i in 1:NA)

        sizes[(lp_total, m_total)] += 1
        push!(vectors[(lp_total, m_total)], copy(v))
    end

    SubsystemBasis{NA}(parent, sort(collect(keys(vectors))), sizes, vectors)
end
