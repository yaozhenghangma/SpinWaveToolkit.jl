struct Bond
    J::Float
    displace::Vector{Float}
end

struct Site
end

struct Lattice
    num_sites::Int
    total_k::Int
end
