function k_space(reciprocal_lattice, num_k)
    total_k = num_k[1] * num_k[2] * num_k[3]
    k = zeros(total_k, 3)
    index = 1
    for ix in range(0, 1, num_k[1]+1)[1:end-1],
        iy in range(0, 1, num_k[2]+1)[1:end-1],
        iz in range(0, 1, num_k[3]+1)[1:end-1]
        k[index, :] = ix * reciprocal_lattice[1] +
                      iy * reciprocal_lattice[2] +
                      iz * reciprocal_lattice[3]
        index += 1
    end
    return k
end

function delta_k_mesh(reciprocal_lattice, num_k)
    total_k = num_k[1] * num_k[2] * num_k[3]
    k = k_space(reciprocal_lattice, num_k)
    ki = permutedims(repeat(k, 1,1, total_k), (1,3,2))
    kj = permutedims(repeat(k, 1,1, total_k), (3,1,2))
    return ki - kj
end
