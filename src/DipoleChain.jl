module DipoleChain

export
    Basis,
    SubsystemBasis,

    sparse_rot,
    sparse_pot,
    dense_x,
    dense_iy,
    dense_z,

    expval1,
    expval2,

    reduced_eigenvalues,
    S_vn,
    S_renyi

include("basis.jl")
include("operators.jl")
include("expectation.jl")
include("entanglement.jl")

end