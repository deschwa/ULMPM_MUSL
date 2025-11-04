# ---------------------------------
# Linear Elastic Isotropic Material
# ---------------------------------

function stress_update!(mp_group::MaterialPointGroup)
    for p_idx in 1:length(mp_group.mass)
        dim = size(mp_group.σ[p_idx], 1)
        E = mp_group.material[p_idx].E
        ν = mp_group.material[p_idx].ν
        λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
        μ = E / (2 * (1 + ν))
        I = I(dim)
        ε_new = 0.5 * (mp_group.L[p_idx] + transpose(mp_group.L[p_idx]))
        tr_ε = tr(ε_new)
        mp_group.σ[p_idx] .+= λ*tr_ε*I + 2*μ*ε_new
    end
end



# --------------------
# Neo-Hookean Material
# --------------------

function stress_update!(mp_group::MaterialPointGroup)
    return
end

