function magnetism(occupation; S=1.0)
    @ein diag_occu[ik, ν] := occupation[ik, ν, ν]
    return S - mean(diag_occu)
end
