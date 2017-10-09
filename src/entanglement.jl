"""
    reduced_eigenvalues{N,NA,NB}(basis::Basis{N}, basis_A::SubsystemBasis{NA}, basis_B::SubsystemBasis{NB}, A_start, wf::AbstractVector{Float64})

Compute the reduced density matrix eigenvalues of state `wf` reduced to
subsystem A.
"""
function reduced_eigenvalues{N,NA,NB}(basis::Basis{N}, basis_A::SubsystemBasis{NA}, basis_B::SubsystemBasis{NB}, A_start, wf::AbstractVector{Float64})
    # Correct subspace sizes.
    N == NA + NB || throw(DomainError())
    # Subsystem is inside the system.
    A_start >= 1 || throw(DomainError())
    A_start + NA - 1 <= N || throw(DomainError())

    v = Array{Int}(2N)
    eigvals = Float64[]

    for block_A in basis_A.blocks
        vectors_A = basis_A.vectors[block_A]

        # Only use blocks that have complementary lp and m.
        block_B = (mod(basis.sym_lp-block_A[1], 2), basis.sym_m-block_A[2])
        vectors_B = basis_B.vectors[block_B]

        mat = zeros(basis_B.sizes[block_B], basis_A.sizes[block_A])

        for col in 1:basis_A.sizes[block_A]
            l_total_col = 0
            for i in 1:2NA
                v[2A_start+i-2] = vectors_A[col][i]
                i % 2 == 1 && (l_total_col += v[2A_start+i-2])
            end

            for row in 1:basis_B.sizes[block_B]
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

                mat[row, col] = wf[basis.lookup[v]]
            end
        end

        append!(eigvals, abs2.(svdvals(mat)))
    end

    abs(sum(eigvals) - 1.0) < 1e-12 || warn("Bad trace: $(abs(sum(eigvals) - 1.0))")

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
