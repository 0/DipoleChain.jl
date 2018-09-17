"""
    reduced_eigenvalues(basis::SingleBlockBasis{N}, basis_A::MultiBlockBasis{NA}, basis_B::MultiBlockBasis{NB}, A_start, wf::AbstractVector{Float64})

Compute the reduced density matrix eigenvalues of state `wf` reduced to
subsystem A.
"""
function reduced_eigenvalues(basis::SingleBlockBasis{N}, basis_A::MultiBlockBasis{NA}, basis_B::MultiBlockBasis{NB}, A_start, wf::AbstractVector{Float64}) where {N,NA,NB}
    N == NA + NB || throw(DomainError((NA, NB), "Incorrect subspace sizes."))
    A_start >= 1 || throw(DomainError((A_start, NA), "Subsystem must be inside system."))
    A_start + NA - 1 <= N || throw(DomainError((A_start, NA), "Subsystem must be inside system."))

    v = Array{Int}(undef, 2N)
    eigvals = Float64[]

    for (label_A, block_A) in basis_A.blocks
        vectors_A = block_A.vectors

        # Only use blocks that have complementary lp and m.
        label_B = BlockLabel(mod(basis.label.lp-label_A.lp, 2), basis.label.m-label_A.m)
        haskey(basis_B.blocks, label_B) || continue

        block_B = basis_B.blocks[label_B]
        vectors_B = block_B.vectors

        mat = zeros(block_B.size, block_A.size)

        for col in 1:block_A.size
            l_total_col = 0
            for i in 1:2NA
                v[2A_start+i-2] = vectors_A[col][i]
                i % 2 == 1 && (l_total_col += v[2A_start+i-2])
            end

            for row in 1:block_B.size
                l_total = l_total_col
                for i in 1:(2A_start-2)
                    v[i] = vectors_B[row][i]
                    i % 2 == 1 && (l_total += v[i])
                end
                for i in (2A_start+2NA-1):2N
                    v[i] = vectors_B[row][i-2NA]
                    i % 2 == 1 && (l_total += v[i])
                end

                l_total <= basis.l_total_max || continue

                mat[row, col] = wf[basis.block.lookup[v]]
            end
        end

        append!(eigvals, abs2.(svdvals(mat)))
    end

    abs(sum(eigvals) - 1.0) < 1e-12 || @warn("Bad trace: $(abs(sum(eigvals) - 1.0))")

    eigvals
end

"""
    S_vn(eigvals::AbstractVector{Float64})

Von Neumann entropy of `eigvals`.
"""
S_vn(eigvals::AbstractVector{Float64}) = -sum(x * log(x) for x in eigvals if x > 0)

"""
    S_renyi(eigvals::AbstractVector{Float64}, alpha=2)

Order-`alpha` Rényi entropy of `eigvals`.
"""
S_renyi(eigvals::AbstractVector{Float64}, alpha=2) = log(sum(eigvals.^alpha)) / (1 - alpha)
