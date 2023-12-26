using Statistics

function magnetism(occupation; S=1.0)
    return S - mean(occupation)
end
