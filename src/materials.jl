# ---------------------------------
# Linear Elastic Isotropic Material
# ---------------------------------

struct LinearElastic <: AbstractMaterial
    E::Float64        # Young's Modulus
    ν::Float64        # Poisson's Ratio
    ρ::Float64        # Density
end

function stress_update!(mp::MaterialPoint{LinearElastic})
    dim = size(mp.σ, 1)
    E = mp.material.E
    ν = mp.material.ν
    λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
    μ = E / (2 * (1 + ν))
    I = I(dim)
    ε_new = 0.5 * (mp.L + transpose(mp.L))
    tr_ε = tr(ε_new)
    mp.σ .+= λ*tr_ε*I + 2*μ*ε_new
    return
end



# --------------------
# Neo-Hookean Material
# --------------------

struct NeoHookean <: AbstractMaterial
    μ::Float64        # Shear Modulus
    λ::Float64       # Lamé's First Parameter
    ρ::Float64        # Density
end

function stress_update!(mp::MaterialPoint{NeoHookean})
    return
end


end