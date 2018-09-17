module DipoleChain

using LinearAlgebra: diag, diagm, I, svdvals, tr
using SparseArrays: sparse

export
    println_result,

    BlockLabel,
    Block,
    AbstractBasis,
    SingleBlockBasis,
    MultiBlockBasis,

    sparse_rot,
    diag_rot,
    sparse_pot,
    dense_pot,
    scale_dense_pot,
    dense_x,
    dense_iy,
    dense_z,

    expval1,
    expval2,

    symmetric_trotter_path,
    wf_pigs,
    energy,
    energy_mixed,

    reduced_eigenvalues,
    S_vn,
    S_renyi

include("util.jl")

include("basis.jl")
include("operators.jl")
include("expectation.jl")
include("path_integrals.jl")
include("entanglement.jl")

end
