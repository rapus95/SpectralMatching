
function IniResponse(a::Vector{Float64},t,m::Float64,tetha::Float64,increment::Float64,
                    Tn::Vector{Float64},wn::Vector{Float64},target::Vector{Float64})

    #Initial set()
    t_index::Vector{Float64}=zeros(length(Tn))
    t_peak::Vector{Float64}=zeros(length(Tn))
    apeak::Vector{Float64}=zeros(length(Tn))
    misfit::Vector{Float64}=zeros(length(Tn))
    u::Vector{Float64}=zeros(length(Tn),length(t))

    #Newmark coefficient
    gamma::Float64=0.5
    betha::Float64=0.25

    @fastmath @inbounds @simd for i=1:length(Tn)
        udot::Float64=0
        @fastmath @inbounds  ks::Float64=m*wn[i]^2
        @fastmath @inbounds c::Float64=2*m*wn[i]*tetha
        @fastmath @inbounds uddot::Float64=(a[1]-ks.*u[i,1])/m

        a1::Float64=1/(betha*increment^2).*m+gamma/(betha*increment).*c
        a2::Float64=1/(betha*increment).*m+(gamma/betha-1).*c
        a3::Float64=(1/2/betha-1).*m+increment*(gamma/2/betha-1).*c
        kh::Float64=ks+a1

        @fastmath @inbounds @simd for j=1:length(t)-1
            @fastmath @inbounds ph::Float64=a[j+1]+a1*u[i,j]+a2*udot+a3*uddot
            @fastmath @inbounds u[i,j+1]=ph/kh
            udoti=udot
            @fastmath @inbounds udot=gamma/(betha*increment).*(u[i,j+1]-u[i,j])+(1-gamma/betha).*udot+increment*(1-gamma/2/betha)*uddot
            @fastmath @inbounds uddot=1/(betha*increment^2).*(u[i,j+1]-u[i,j])-1/(betha*increment).*udoti-(1/2/betha-1)*uddot
        end
        
        @fastmath @inbounds abslt::Vector{Float64}=abs.(u[i,:]).*wn[i].^2;              #Absolute pseudo acceleration response [psa]
        R::Float64=maximum(abslt);                           #Max of absolute psa
        @fastmath @inbounds t_index[i]=findall(x->x==R,abslt)[1];              #Index of max absolute psa on time series
        @fastmath @inbounds t_peak[i]=(t_index[i]-1).*increment;    #Time of peak psa
        @fastmath @inbounds apeak[i]=u[i,convert(Int64,t_index[i])].*wn[i]^2;      #True value of max peak
    end

    #Calculate misfit
    @fastmath @inbounds @simd for i=1:length(Tn)
        if apeak[i]<0
            @fastmath @inbounds target[i]=-abs(target[i]);
        end   
        @fastmath @inbounds misfit[i]=target[i]-apeak[i]; 
    end


    return u,t_peak,t_index,apeak,misfit,target;
end