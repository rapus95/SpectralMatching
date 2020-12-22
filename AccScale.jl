function AccScale(apeak::Vector{Float64},target::Vector{Float64},acc::Vector{Float64})
    scale_factor::Float64 = sum(target .*apeak)/sum(apeak.^2)
    return scale_factor*acc
end
