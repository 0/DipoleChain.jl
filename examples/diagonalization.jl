#!/usr/bin/env julia

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "../src"))
using DipoleChain

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "-R"
        metavar = "R"
        help = "separation distance"
        arg_type = Float64
        required = true
    "-N"
        metavar = "N"
        help = "number of rotors"
        arg_type = Int
        required = true
    "--l-max"
        metavar = "L"
        help = "local basis truncation"
        arg_type = Int
        required = true
    "--l-total_max"
        metavar = "L"
        help = "many-body basis truncation (default: no truncation)"
        arg_type = Int
    "--sym-lp"
        metavar = "LP"
        help = "total l parity of symmetry block"
        arg_type = Int
        default = 0
    "--sym-m"
        metavar = "M"
        help = "total m value of symmetry block"
        arg_type = Int
        default = 0
    "--A-start"
        metavar = "jA"
        help = "first rotor of subsystem A"
        arg_type = Int
        default = 1
    "--A-size"
        metavar = "NA"
        help = "length of subsystem A"
        arg_type = Int
        default = 1
end
c = parse_args(ARGS, s, as_symbols=true)

R = c[:R]
N = c[:N]
l_max = c[:l_max]
l_total_max = c[:l_total_max]
if l_total_max === nothing
    l_total_max = l_max * N
end
sym_lp = c[:sym_lp]
sym_m = c[:sym_m]
A_start = c[:A_start]
A_size = c[:A_size]


println("[ ] Constructing basis.")
@time basis = SingleBlockBasis(N, l_max, l_total_max, sym_lp, sym_m)
println("[+] Constructed basis: $(basis.size).")

println("[ ] Constructing matrix elements.")
@time K = sparse_rot(basis)
@time V = sparse_pot(basis)
@time H = K + V/R^3
# Number of non-zero elements.
println("[+] Constructed matrix elements: $(nnz(H)).")

println("[ ] Diagonalizing.")
@time d = eigs(H, which=:SR, nev=1)
E0 = d[1][1]
wf = d[2][:,1]
# Number of iterations.
println("[+] Done diagonalizing: $(d[4]).")

# Ground state energy.
println_result("E0 = $(E0)")


println("[ ] Calculating correlations.")
@time for (op_name, op_f) in [("x", dense_x), ("z", dense_z)]
    println_result("<$(op_name)_i $(op_name)_j>")
    op = op_f(basis)
    for i in 1:N
        for j in 1:(i-1)
            @printf("%18s ", "")
        end
        @printf(" % .15f", expval1(basis, op*op, i, wf))
        for j in (i+1):N
            @printf(" % .15f", expval2(basis, op, i, op, j, wf))
        end
        println()
    end
end
println("[+] Calculated correlations.")


println("[ ] Constructing subsystem bases.")
@time basis_A = MultiBlockBasis(A_size, l_max, l_total_max)
@time basis_B = MultiBlockBasis(N-A_size, l_max, l_total_max)
println("[+] Constructed subsystem bases: $(basis_A.size), $(basis_B.size).")

println("[ ] Calculating reduced eigenvalues.")
@time eigvals = reduced_eigenvalues(basis, basis_A, basis_B, A_start, wf)
println("[+] Calculated reduced eigenvalues.")

# Von Neumann entropy.
println_result("SvN = $(S_vn(eigvals))")
# Order-2 Rényi entropy.
println_result("S2 = $(S_renyi(eigvals, 2))")
