"""
    eye(n::Int)

Create an `n` by `n` identity matrix.
"""
eye(n::Int) = Matrix{Float64}(I, n, n)

"""
    eye(n1::Int, n2::Int)

Create an `n1` by `n2` identity matrix.
"""
eye(n1::Int, n2::Int) = Matrix{Float64}(I, n1, n2)


"""
    println_result(msg...)

Output a result so that it stands out.
"""
println_result(msg...) = printstyled(msg..., "\n"; color=:green, bold=true)
