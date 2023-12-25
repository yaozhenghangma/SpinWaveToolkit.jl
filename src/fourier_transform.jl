using OMEinsum

function k_transform(k, bonds)
    Jk = zeros(size(k)[1])
    for bond in bonds
        @ein tmp[i] := k[i, k] * bond.displace[k]
        Jk += tmp * bond.J
    end
    return Jk
end

function k_mesh_transform(k_mesh, bonds)
    Jkk = zeros(size(k)[1], size(k)[2])
    for bond in bonds
        @ein tmp[i, j] := k_mesh[i, j, k] * bond.displace[k]
    end
    return Jkk
end
