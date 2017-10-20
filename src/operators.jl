function sparse_rot{N}(block::Block{N})
    rows = Int[]
    elems = Float64[]

    for row in 1:block.size
        elem = 0.0

        for i in 1:N
            l = block.vectors[row][2i-1]
            elem += l*(l+1)
        end

        push!(rows, row)
        push!(elems, elem)
    end

    sparse(rows, rows, elems, block.size, block.size)
end

"""
    sparse_rot{N}(basis::SingleBlockBasis{N})

Compute a sparse representation of the dimensionless rotational energy matrix.
"""
sparse_rot{N}(basis::SingleBlockBasis{N}) = sparse_rot(basis.block)

"""
    diag_rot{N}(basis::SingleBlockBasis{N})

Compute a diagonal representation of the dimensionless rotational energy
matrix.
"""
diag_rot{N}(basis::SingleBlockBasis{N}) = diag(sparse_rot(basis.block))

"""
    diag_rot{N}(basis::MultiBlockBasis{N})

Compute a diagonal representation of the dimensionless rotational energy
matrix.
"""
function diag_rot{N}(basis::MultiBlockBasis{N})
    result = Dict{BlockLabel,Vector{Float64}}()

    for (label, block) in basis.blocks
        result[label] = diag(sparse_rot(block))
    end

    result
end


function C(l::Int, lp::Int, m::Int)
    if lp == l-1
        sqrt((l-m)*(l+m)/((2l-1)*(2l+1)))
    elseif lp == l+1
        sqrt((lp-m)*(lp+m)/((2lp-1)*(2lp+1)))
    end
end

function Sp(l::Int, lp::Int, m::Int)
    if lp == l-1
        -sqrt((l-m-1)*(l-m)/((2l-1)*(2l+1)))
    elseif lp == l+1
        sqrt((lp+m)*(lp+m+1)/((2lp-1)*(2lp+1)))
    end
end

function Sm(l::Int, lp::Int, m::Int)
    if lp == l-1
        sqrt((l+m-1)*(l+m)/((2l-1)*(2l+1)))
    elseif lp == l+1
        -sqrt((lp-m)*(lp-m+1)/((2lp-1)*(2lp+1)))
    end
end

function pot_elem(l1::Int, m1::Int, l2::Int, m2::Int, l1p::Int, m1p::Int, l2p::Int, m2p::Int)
    if m1p == m1-1 && m2p == m2+1
        Sm(l1, l1p, m1) * Sp(l2, l2p, m2)
    elseif m1p == m1+1 && m2p == m2-1
        Sp(l1, l1p, m1) * Sm(l2, l2p, m2)
    elseif m1p == m1 && m2p == m2
        -4.0 * C(l1, l1p, m1) * C(l2, l2p, m2)
    end
end


function sparse_pot{N}(basis::AbstractBasis{N}, block::Block{N})
    rows = Int[]
    cols = Int[]
    elems = Float64[]

    w = Array{Int}(2N)

    for col in 1:block.size
        v = block.vectors[col]
        l_total = 0
        for i in 1:N
            l_total += v[2i-1]
        end

        # All rotor pairs.
        for i in 1:(N-1)
            for j in (i+1):N
                # All potentially allowed transitions.
                for delta_li in [-1, 1]
                    for delta_lj in [-1, 1]
                        for delta_m in -1:1
                            w .= v
                            w[2i-1] += delta_li
                            w[2i] += delta_m
                            w[2j-1] += delta_lj
                            w[2j] -= delta_m

                            # Bounds on l.
                            0 <= w[2i-1] <= basis.l_max || continue
                            0 <= w[2j-1] <= basis.l_max || continue
                            l_total + delta_li + delta_lj <= basis.l_total_max || continue
                            # Bounds on m.
                            -w[2i-1] <= w[2i] <= w[2i-1] || continue
                            -w[2j-1] <= w[2j] <= w[2j-1] || continue

                            row = block.lookup[w]
                            elem = 0.5 * pot_elem(v[2i-1], v[2i], v[2j-1], v[2j], w[2i-1], w[2i], w[2j-1], w[2j]) / (j-i)^3

                            push!(rows, row)
                            push!(cols, col)
                            push!(elems, elem)
                        end
                    end
                end
            end
        end
    end

    sparse(rows, cols, elems, block.size, block.size)
end

"""
    sparse_pot{N}(basis::SingleBlockBasis{N})

Compute a sparse representation of the dimensionless potential energy matrix
for dipole-dipole interactions between all rotor pairs.
"""
sparse_pot{N}(basis::SingleBlockBasis{N}) = sparse_pot(basis, basis.block)

"""
    dense_pot{N}(basis::SingleBlockBasis{N})

Compute a dense representation of the dimensionless potential energy matrix for
dipole-dipole interactions between all rotor pairs.
"""
dense_pot{N}(basis::SingleBlockBasis{N}) = full(sparse_pot(basis, basis.block))

"""
    dense_pot{N}(basis::MultiBlockBasis{N})

Compute a dense representation of the dimensionless potential energy matrix for
dipole-dipole interactions between all rotor pairs.
"""
function dense_pot{N}(basis::MultiBlockBasis{N})
    result = Dict{BlockLabel,Matrix{Float64}}()

    for (label, block) in basis.blocks
        result[label] = full(sparse_pot(basis, block))
    end

    result
end


"""
    scale_dense_pot(pot::Matrix{Float64}, k::Float64)

Scale the matrix elements of `pot` in place by a factor of `k`.
"""
function scale_dense_pot(pot::Matrix{Float64}, k::Float64)
    pot .*= k
end

"""
    scale_dense_pot(pot::Dict{BlockLabel,Matrix{Float64}}, k::Float64)

Scale the matrix elements of `pot` in place by a factor of `k`.
"""
function scale_dense_pot(pot::Dict{BlockLabel,Matrix{Float64}}, k::Float64)
    for block in values(pot)
        block .*= k
    end
end


"""
    dense_x{N}(basis::SingleBlockBasis{N})

Compute a dense representation of the single-rotor x operator.
"""
function dense_x{N}(basis::SingleBlockBasis{N})
    l_max = basis.l_max

    result = zeros(Float64, (l_max+1)^2, (l_max+1)^2)

    for l in 0:l_max
        for m in (-l):l
            col = l*(l+1) + m + 1

            for lp in 0:l_max
                for mp in (-lp):lp
                    row = lp*(lp+1) + mp + 1

                    if lp == l+1 && mp == m+1
                        result[row, col] -= 0.5 * sqrt((l+m+1)*(l+m+2)/((2l+1)*(2l+3)))
                    elseif lp == l-1 && mp == m+1
                        result[row, col] += 0.5 * sqrt((l-m-1)*(l-m)/((2l-1)*(2l+1)))
                    elseif lp == l+1 && mp == m-1
                        result[row, col] += 0.5 * sqrt((l-m+1)*(l-m+2)/((2l+1)*(2l+3)))
                    elseif lp == l-1 && mp == m-1
                        result[row, col] -= 0.5 * sqrt((l+m-1)*(l+m)/((2l-1)*(2l+1)))
                    end
                end
            end
        end
    end

    result
end

"""
    dense_iy{N}(basis::SingleBlockBasis{N})

Compute a dense representation of the single-rotor iy operator.
"""
function dense_iy{N}(basis::SingleBlockBasis{N})
    l_max = basis.l_max

    result = zeros(Float64, (l_max+1)^2, (l_max+1)^2)

    for l in 0:l_max
        for m in (-l):l
            col = l*(l+1) + m + 1

            for lp in 0:l_max
                for mp in (-lp):lp
                    row = lp*(lp+1) + mp + 1

                    if lp == l+1 && mp == m+1
                        result[row, col] -= 0.5 * sqrt((l+m+1)*(l+m+2)/((2l+1)*(2l+3)))
                    elseif lp == l-1 && mp == m+1
                        result[row, col] += 0.5 * sqrt((l-m-1)*(l-m)/((2l-1)*(2l+1)))
                    elseif lp == l+1 && mp == m-1
                        result[row, col] -= 0.5 * sqrt((l-m+1)*(l-m+2)/((2l+1)*(2l+3)))
                    elseif lp == l-1 && mp == m-1
                        result[row, col] += 0.5 * sqrt((l+m-1)*(l+m)/((2l-1)*(2l+1)))
                    end
                end
            end
        end
    end

    result
end

"""
    dense_z{N}(basis::SingleBlockBasis{N})

Compute a dense representation of the single-rotor z operator.
"""
function dense_z{N}(basis::SingleBlockBasis{N})
    l_max = basis.l_max

    result = zeros(Float64, (l_max+1)^2, (l_max+1)^2)

    for l in 0:l_max
        for m in (-l):l
            col = l*(l+1) + m + 1

            for lp in 0:l_max
                for mp in (-lp):lp
                    row = lp*(lp+1) + mp + 1

                    if lp == l+1 && mp == m
                        result[row, col] += sqrt((l-m+1)*(l+m+1)/((2l+1)*(2l+3)))
                    elseif lp == l-1 && mp == m
                        result[row, col] += sqrt((l-m)*(l+m)/((2l-1)*(2l+1)))
                    end
                end
            end
        end
    end

    result
end
