function Zeropad(zeropad::Int64,arec::Array{Float64},dt::Float64)
    pad::Int64 = convert(Int64,zeropad/dt);
    apad::Array{Float64} =zeros(convert(Int64,2*pad+length(arec)));
    apad[pad+1: pad+length(arec)] = arec;
    t::Array{Float64}=(0:length(apad)-1).*dt;
    return apad',t'
end