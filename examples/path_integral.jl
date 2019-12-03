#!/usr/bin/env julia

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
    "--beta"
        metavar = "B"
        help = "reciprocal temperature"
        arg_type = Float64
    "--tau"
        metavar = "T"
        help = "reciprocal temperature per link"
        arg_type = Float64
    "-P"
        metavar = "P"
        help = "number of beads"
        arg_type = Int
        required = true
    "--pigs"
        help = "use PIGS"
        action = :store_true
    "--A-start"
        metavar = "jA"
        help = "first rotor of subsystem A (PIGS only)"
        arg_type = Int
        default = 1
    "--A-size"
        metavar = "NA"
        help = "length of subsystem A (PIGS only)"
        arg_type = Int
        default = 1
end
c = parse_args(ARGS, s, as_symbols=true)

R = c[:R]
N = c[:N]
l_max = c[:l_max]
l_total_max = c[:l_total_max]
if isnothing(l_total_max)
    l_total_max = l_max * N
end
beta = c[:beta]
tau = c[:tau]
if isnothing(beta) && isnothing(tau) || !isnothing(beta) && !isnothing(tau)
    println("Exactly one of --beta or --tau must be set.")
    exit()
end
P = c[:P]
pigs = c[:pigs]
A_start = c[:A_start]
A_size = c[:A_size]

if pigs
    if P % 2 == 0
        println("P must be odd for PIGS.")
        exit()
    end
    num_links = P - 1
else
    num_links = P
end

if isnothing(beta)
    beta = tau * num_links
elseif isnothing(tau)
    tau = beta / num_links
end


println("[ ] Constructing basis.")
@time if pigs
    # Always use the ground state symmetry block for PIGS.
    sym_lp = 0
    sym_m = 0
    basis = SingleBlockBasis(N, l_max, l_total_max, sym_lp, sym_m)
else
    basis = MultiBlockBasis(N, l_max, l_total_max)
end
println("[+] Constructed basis: $(basis.size).")

println("[ ] Constructing matrix elements.")
@time K = diag_rot(basis)
@time V = dense_pot(basis)
@time scale_dense_pot(V, 1.0/R^3)
println("[+] Constructed matrix elements.")

println("[ ] Constructing path.")
@time if pigs
    # Trial function is free rotor ground state.
    trial_idx = basis.block.lookup[[0 for _ in 1:2N]]

    half_path = symmetric_trotter_path(div(num_links, 2), tau, K, V)
    path = half_path * half_path
    wf = wf_pigs(half_path, trial_idx)
else
    path = symmetric_trotter_path(num_links, tau, K, V)
end
println("[+] Constructed path.")

println("[ ] Calculating energies.")
@time if pigs
    E0, (avgK, avgV) = energy_mixed(path, trial_idx, K, V)
    # For a spherically symmetric trial function, the kinetic contribution to
    # the mixed estimator should vanish.
    avgK == 0.0 || @warn("avgK == $(avgK)")
else
    Z, E, (avgK, avgV) = energy(path, K, V)
end
println("[+] Calculated energies.")

if pigs
    # Ground state energy.
    println_result("E0 = $(E0)")
else
    # Partition function.
    println_result("Z = $(Z)")
    # Rotational energy.
    println_result("K = $(avgK)")
    # Potential energy.
    println_result("V = $(avgV)")
    # Total energy.
    println_result("E = $(E)")
end


if pigs
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
end
