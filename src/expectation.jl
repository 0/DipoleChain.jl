"""
    expval1(basis::SingleBlockBasis, op::AbstractMatrix, i::Int, wf::AbstractVector{Float64})

Calculate the expectation value of the one-body operator `op` at site `i` in
the state `wf`.
"""
function expval1(basis::SingleBlockBasis{N}, op::AbstractMatrix, i::Int, wf::AbstractVector{Float64}) where {N}
    w = Array{Int}(undef, 2N)

    result = 0.0

    for col in 1:basis.block.size
        v = basis.block.vectors[col]

        li = v[2i-1]
        mi = v[2i]
        coli = li*(li+1) + mi + 1

        rowi = 0
        for lip in 0:basis.l_max
            for mip in (-lip):lip
                rowi += 1
                op[rowi, coli] == 0.0 && continue

                w .= v
                w[2i-1] = lip
                w[2i] = mip

                haskey(basis.block.lookup, w) || continue

                result += wf[basis.block.lookup[w]] * wf[col] * op[rowi, coli]
            end
        end
    end

    result
end

"""
    expval2(basis::SingleBlockBasis, op_i::AbstractMatrix, i::Int, op_j::AbstractMatrix, j::Int, wf::AbstractVector{Float64})

Calculate the expectation value of the two-body operator composed of `op_i` at
site `i` and `op_j` at site `j` in the state `wf`.
"""
function expval2(basis::SingleBlockBasis{N}, op_i::AbstractMatrix, i::Int, op_j::AbstractMatrix, j::Int, wf::AbstractVector{Float64}) where {N}
    w = Array{Int}(undef, 2N)

    result = 0.0

    for col in 1:basis.block.size
        v = basis.block.vectors[col]

        li = v[2i-1]
        mi = v[2i]
        coli = li*(li+1) + mi + 1
        lj = v[2j-1]
        mj = v[2j]
        colj = lj*(lj+1) + mj + 1

        rowi = 0
        for lip in 0:basis.l_max
            for mip in (-lip):lip
                rowi += 1
                op_i[rowi, coli] == 0.0 && continue

                rowj = 0
                for ljp in 0:basis.l_max
                    for mjp in (-ljp):ljp
                        rowj += 1
                        op_j[rowj, colj] == 0.0 && continue

                        w .= v
                        w[2i-1] = lip
                        w[2i] = mip
                        w[2j-1] = ljp
                        w[2j] = mjp

                        haskey(basis.block.lookup, w) || continue

                        result += wf[basis.block.lookup[w]] * wf[col] * op_i[rowi, coli] * op_j[rowj, colj]
                    end
                end
            end
        end
    end

    result
end
