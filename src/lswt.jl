using LinearAlgebra

KB = 8.6173e-2
gs = 2.02
μB = 0.057884
ħ = 1.0

function lswt(lattice, T; S=1.0, A=0, B=0)
    kbt = 1.0/(T*KB)
    @. H = S*lattice.Jk
    for i in 1:lattice.num_sites
        Jki = sum(lattice.Jk[1, i, :])
        @. H[:, i, i] += -S*Jki -2*S*A + gs*μB*B
    end

    energies = zeros(lattice.total_k, lattice.num_sites)
    for ik in 1:lattice.total_k
        energies[ik, :] = eigvals(H[ik, :, :])
    end
    @. occupation = 1.0/(exp(energies*ħ*kbt) -1.0)
    occupation = sum(occupation, dims=2)
    return energies, occupation
end
