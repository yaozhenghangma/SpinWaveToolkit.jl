struct Bond
    J::Float64
    displace::Vector{Float64}
    linked_site::Int
end

mutable struct Site
    index::Int64
    bonds::Vector{Bond}
end

struct Lattice
    num_sites::Int64
    total_k::Int64
    Jk::Array{Float64, 3}
    Jkk::Array{Float64, 4}
end
