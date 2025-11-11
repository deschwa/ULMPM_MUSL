using LinearAlgebra
# ---------------------------------
# Linear Elastic Isotropic Material
# ---------------------------------

function stress_update!(mp_group::MaterialPointGroup, dt::Float64)
    dim = size(mp_group.σ[1], 1)
    E = mp_group.material.E
    ν = mp_group.material.ν
    λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
    μ = E / (2 * (1 + ν))
    I_dim = Matrix(I(dim))
    
    for p_idx in 1:length(mp_group.mass)
        ε_new = 0.5 * dt * (mp_group.L[p_idx] + transpose(mp_group.L[p_idx]))
        tr_ε = tr(ε_new)
        mp_group.σ[p_idx] .+= λ*tr_ε*I_dim + 2*μ*ε_new
    end
end



# --------------------
# Neo-Hookean Material
# --------------------
