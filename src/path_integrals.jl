function expmk(k::Float64, op::AbstractVector{Float64})
    diagm(0 => exp.(k .* op))
end

function expmk(k::Float64, op::AbstractMatrix{Float64})
    exp(k .* op)
end


"""
    symmetric_trotter_path(num_links::Int, tau::Float64, ops::AbstractArray...)

Compute an approximation of `exp(-num_links*tau*sum(ops))` using a symmetric
Trotter decomposition.

The first term in `ops` is placed in the middle of the decomposition.
"""
function symmetric_trotter_path(num_links::Int, tau::Float64, ops::AbstractArray...)
    step = eye(size(ops[1])...)

    for i in length(ops):-1:2
        step *= expmk(-tau/2, ops[i])
    end
    step *= expmk(-tau, ops[1])
    for i in 2:length(ops)
        step *= expmk(-tau/2, ops[i])
    end

    step^num_links
end

"""
    symmetric_trotter_path(num_links::Int, tau::Float64, ops::Dict...)

Compute an approximation of `exp(-num_links*tau*sum(ops))` using a symmetric
Trotter decomposition.

The first term in `ops` is placed in the middle of the decomposition.
"""
function symmetric_trotter_path(num_links::Int, tau::Float64, ops::Dict...)
    path = Dict{BlockLabel,Matrix{Float64}}()

    for label in keys(ops[1])
        step = eye(size(ops[1][label])...)

        for i in length(ops):-1:2
            step *= expmk(-tau/2, ops[i][label])
        end
        step *= expmk(-tau, ops[1][label])
        for i in 2:length(ops)
            step *= expmk(-tau/2, ops[i][label])
        end

        path[label] = step^num_links
    end

    path
end


"""
    wf_pigs(half_path::Matrix{Float64}, trial_idx::Int)

Extract the propagated ground state wavefunction from the middle of a PIGS
path.

To access the middle of the path, we require only one half of the path in
`half_path`. The half path is multiplied by a trial function, which is the
basis vector specified by `trial_idx`.
"""
function wf_pigs(half_path::Matrix{Float64}, trial_idx::Int)
    wf = half_path[:, trial_idx]
    wf ./= sqrt(sum(wf.^2))

    wf
end


mult_by_op(op1::AbstractMatrix{Float64}, op2::AbstractVector{Float64}) = op1 * diagm(0 => op2)
mult_by_op(op1::AbstractMatrix{Float64}, op2::AbstractMatrix{Float64}) = op1 * op2


"""
    energy(path::Dict, ops::Dict...)

Compute the thermal energy of a finite temperature `path` for a Hamiltonian
with terms `ops`.
"""
function energy(path::Dict, ops::Dict...)
    Z = 0.0
    avgs = zeros(Float64, length(ops))

    for label in keys(path)
        Z += tr(path[label])

        for i in 1:length(ops)
            avgs[i] += tr(mult_by_op(path[label], ops[i][label]))
        end
    end

    avgs /= Z
    E = sum(avgs)

    Z, E, avgs
end


"""
    energy_mixed(path::Matrix{Float64}, trial_idx::Int, ops::AbstractArray...)

Compute the ground state energy of a PIGS `path` using the mixed estimator for
a Hamiltonian with terms `ops`.

The path is multiplied on both ends by a trial function, which is the basis
vector specified by `trial_idx`.
"""
function energy_mixed(path::Matrix{Float64}, trial_idx::Int, ops::AbstractArray...)
    avgs = Array{Float64}(undef, length(ops))

    for i in 1:length(ops)
        avgs[i] = mult_by_op(path, ops[i])[trial_idx, trial_idx]
    end

    avgs /= path[trial_idx, trial_idx]
    E0 = sum(avgs)

    E0, avgs
end
