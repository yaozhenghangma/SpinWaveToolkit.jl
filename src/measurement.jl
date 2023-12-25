using Statistics

function magnetism(occupation, S)
    return S - mean(occupation)
end
