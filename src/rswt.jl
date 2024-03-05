#function rswt(lattice, T; S=1.0, A=0, B=0, δ=1e-6, max_loop=1e4, λ=1.0)
#    energies, occupation = lswt(lattice, T, S=S, A=A, B=B)
#
#    kbt = 1.0/(T*KB)
#    loop_num = 0
#    Δ = 1.0
#    while Δ > δ && loop_num < max_loop
#        loop_num += 1
#        average_occupation = mean(occupation) / lattice.num_sites
#        H = @. lattice.Jk * (S - average_occupation)
#        for i in 1:lattice.num_sites
#            Jki = sum(lattice.Jk[1, i, :])
#            @. H[:, i, i] += -S*Jki - 2*S*A + gs*μB*B + 4*A*average_occupation + A
#            @ein tmp1[k2] := lattice.Jkk[i, k1, k2] * occupation[k1, 1]
#            tmp2 = mean(lattice.Jk[:, i, i] .* occupation)
#            @. H[:, i, i] += tmp1/lattice.total_k + tmp2 + Jki * average_occupation
#        end
#
#        for ik in 1:lattice.total_k
#            energies[ik, :] = abs.(eigvals(H[ik, :, :]))*λ + energies[ik, :]*(1-λ)
#        end
#        new_occupation = @. 1.0/(exp(energies*ħ*kbt) -1.0)
#        new_occupation = sum(new_occupation, dims=2)
#        Δ = maximum(abs.(occupation-new_occupation))
#        occupation = deepcopy(new_occupation)
#        loop_num += 1
#    end
#    return energies, occupation, loop_num
#end

function rswt(lattice, T; S=1.0, A=0, B=0, δ=1e-6, max_loop=1e4, λ=1.0)
    energies, occupation = lswt(lattice, T, S=S, A=A, B=B)
    occupation = zeros(size(occupation))

    kbt = 1.0/(T*KB)
    loop_num = 0
    Δ = 1.0
    while Δ > δ && loop_num < max_loop
        loop_num += 1
        average_occupation = mean(occupation) / lattice.num_sites

        H = @. lattice.Jk * S

        sum_occu = diag(dropdims(sum(occupation, dims=1), dims=1))
        sum_occu_matrix = permutedims(repeat(sum_occu, 1, lattice.num_sites), (1,2))
        sum_occu_matrix += permutedims(repeat(sum_occu, 1, lattice.num_sites), (2,1))
        sum_occu_matrix .*= 1.0/(2*lattice.total_k)
        sum_occu_matrix = permutedims(repeat(sum_occu_matrix, 1,1,lattice.total_k), (3,1,2))
        H -= sum_occu_matrix .* lattice.Jk

        occu_matrix = permutedims(repeat(occupation, 1,1,1,lattice.total_k), (4,1,2,3))
        tmp = occu_matrix .* lattice.Jkk
        tmp = dropdims(sum(tmp, dims=2), dims=2)
        H += tmp .* (1.0/lattice.total_k)

        # diagonal term
        H_diag = fill(gs*μB*B + (1-2*S)*A, lattice.num_sites)
        H_diag += sum_occu .* (4*A/lattice.total_k)
        Jk0 = lattice.Jk[1, :, :]
        H_diag -= dropdims(sum(Jk0, dims=2), dims=2) .* S
        @ein temp1[i] := sum_occu[j] * Jk0[i, j]
        H_diag += temp1 .* (1.0/lattice.total_k)
        @ein temp2[i, j] := lattice.Jnk[ik, i, l] * occupation[ik, l, j]
        @ein temp3[i, j] := lattice.Jnk[ik, l, i] * occupation[ik, j, l]
        H_diag -= (diag(temp2) + diag(temp3)) .* (0.5/lattice.total_k)
        H_diag_matrix = diagm(H_diag)
        H += permutedims(repeat(H_diag_matrix, 1,1,lattice.total_k), (3,1,2))



        #for i in 1:lattice.num_sites
        #    Jki = sum(lattice.Jk[1, i, :])
        #    @. H[:, i, i] += -S*Jki - 2*S*A + gs*μB*B + 4*A*average_occupation + A
        #    @ein tmp1[k2] := lattice.Jkk[i, k1, k2] * occupation[k1, 1]
        #    tmp2 = mean(lattice.Jk[:, i, i] .* occupation)
        #    @. H[:, i, i] += tmp1/lattice.total_k + tmp2 + Jki * average_occupation
        #end

        num_branch = lattice.num_sites
        Λ = zeros(lattice.total_k, num_branch, lattice.num_sites)
        for ik in 1:lattice.total_k
            energies[ik, :], Λ[ik, :, :] = eigen(H[ik, :, :])
        end
        energies = abs.(energies)
        fermi_func = @. 1.0/(exp(energies*ħ*kbt) -1.0)
        @ein new_occu[ik, ν1, ν2] := conj(Λ)[ik, ν1, σ] * Λ[ik, ν2, σ] * fermi_func[ik, σ]
        Δ = maximum(abs.(occupation-new_occu))
        occupation = deepcopy(new_occu)
        loop_num += 1
    end
    return energies, occupation, loop_num
end
