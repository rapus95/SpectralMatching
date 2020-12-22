function AccScale(apeak::Array{Float64},target::Array{Float64},acc::Array{Float64})
    scale_factor::Float64 = sum(target .*apeak)/sum(apeak.^2)
    return scale_factor*acc
end
