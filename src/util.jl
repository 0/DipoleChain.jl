"""
    println_result(msg...)

Output a result so that it stands out.
"""
println_result(msg...) = printstyled(msg..., "\n"; color=:green, bold=true)
