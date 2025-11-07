# using Threads
using LinearAlgebra


"""
Compute the linear shape functions and their gradients for a given relative position.
"""
function shape_function(r_rel::MVector{2, Float64}, dx::Float64, dy::Float64)
           
    N_x = max(1 - abs(r_rel[1]) / dx, 0.0)
    N_y = max(1 - abs(r_rel[2]) / dy, 0.0)
    N_I = N_x * N_y


    dN_xdx = - sign(r_rel[1]) / dx
    dN_ydy = - sign(r_rel[2]) / dy
    ∇N_I = [dN_xdx * N_y,
            dN_ydy * N_x]


    return (N_I, ∇N_I)
end