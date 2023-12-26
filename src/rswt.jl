using Statistics
using OMEinsum

include("lswt.jl")

function rswt(lattice, T; S=1.0, A=0, B=0, δ=1e-6, max_loop=1e4, λ=1.0)
    energies, occupation = lswt(lattice, T, S=S, A=A, B=B)

    kbt = 1.0/(T*KB)
    loop_num = 0
    Δ = 1.0
    while Δ > δ && loop_num < max_loop
        loop_num += 1
        average_occupation = mean(occupation) / lattice.num_sites
        H = @. lattice.Jk * (S - average_occupation)
        for i in 1:lattice.num_sites
            Jki = sum(lattice.Jk[1, i, :])
            @. H[:, i, i] += -S*Jki - 2*S*A + gs*μB*B + 4*A*average_occupation + A
            @ein tmp1[k2] := lattice.Jkk[i, k1, k2] * occupation[k1, 1]
            tmp2 = mean(lattice.Jk[:, i, i] .* occupation)
            @. H[:, i, i] += tmp1/lattice.total_k + tmp2 + Jki * average_occupation
        end

        for ik in 1:lattice.total_k
            energies[ik, :] = eigvals(H[ik, :, :])*λ + energies[ik, :]*(1-λ)
        end
        new_occupation = @. 1.0/(exp(energies*ħ*kbt) -1.0)
        new_occupation = sum(new_occupation, dims=2)
        Δ = maximum(abs.(occupation-new_occupation))
        occupation = deepcopy(new_occupation)
        loop_num += 1
    end
    return energies, occupation, loop_num
end
