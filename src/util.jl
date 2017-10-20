"""
    println_result(msg...)

Output a result so that it stands out.
"""
println_result(msg...) = print_with_color(:green, msg..., "\n", bold=true)
