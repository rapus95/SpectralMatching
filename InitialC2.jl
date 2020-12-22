function InitialC2(T::Vector{Float64},w::Vector{Float64},m::Float64,tetha::Float64,increment::Float64,
                    t_peak::Vector{Float64},t_index::Vector{Float64},t,Cgain::Float64,
                    target::Vector{Float64},misfit::Vector{Float64},amod::Vector{Float64},CoffDiag::Float64,wavmag::Float64)
    #This routine calculate the C matrix & initial b vector

    #initial set()
    b::Vector{Float64}=ones(length(T))
    tw_index::Vector{Int64}=zeros(length(T))
    tw_peak::Vector{Float64}=zeros(length(T))
    dTj::Vector{Float64}=zeros(length(T));  
    wav::Vector{Float64}=zeros(length(T),length(t))
    gamma::Float64=0.5
    betha::Float64=0.25
    uw::Vector{Float64}=Vector{Float64}(undef,length(t))
    @. wj::Vector{Float64}=w*sqrt(1-tetha^2);  
    @. tj::Vector{Float64}=t-t_peak+dTj;  
    @. gf::Vector{Float64}=1.178*(1/T)^-0.93
    @. ks::Vector{Int64}=m*w^2
    @. c::Vector{Int64}=2*m*w*tetha
    @. a1::Vector{Int64}=1/(betha*increment^2)*m+gamma/(betha*increment)*c
    @. a2::Vector{Int64}=1/(betha*increment)*m+(gamma/betha-1)*c
    @. a3::Vector{Int64}=(1/2/betha-1)*m+increment*(gamma/2/betha-1)*c
    @. kh::Vector{Int64}=ks+a1
    for i=1:length(T) #structure
        fill!(0,uw)
        uwdot::Float64=0
        @. wav[i,:]=wavmag*cos(wj[i]*tj[i])*exp(-(tj[i]/gf[i])^2)*m*b[i]
        
        @. uwddot::Float64=(wav[i,1]-c[i]*uwdot-ks[i]*uw[1])/m        
        for j=1:length(t)-1
            ph::Float64=wav[i,j+1]+a1[i]*uw[j]+a2[i]*uwdot+a3[i]*uwddot
            uw[j+1]=ph/kh[i]
            uwdoti::Float64=uwdot
            uwdot=gamma/(betha*increment).*(uw[j+1]-uw[j])+(1-gamma/betha).*uwdot+increment*(1-gamma/2/betha)*uwddot
            uwddot=1/(betha*increment^2).*(uw[j+1]-uw[j])-1/(betha*increment).*uwdoti-(1/2/betha-1)*uwddot
        end
        # Find peak time of each wavelet   
        abslt::Vector{Float64}=abs.(uw).*w[i].^2
        R::Float64=maximum(abslt)

        tw_index[i]=findall(x->x==R,abslt)[1];
        tw_peak[i]=(tw_index[i]-1).*increment; 

        dTj[i]=tw_peak[i]-t_peak[i]
        @. tj=t-t_peak[i]+dTj[i]
        wav[i,:]=wavmag.*cos.(wj[i].*tj).*exp.(-(tj./gf[i]).^2).*m.*b[i]
    end
   
    tw_peak=zeros(length(T),length(T))
    misfitini::Vector{Float64}=zeros(length(T))
    C::Vector{Float64}=ones(length(T),length(T))
    wavtot::Vector{Float64}=zeros(length(t))

    for i=1:length(T) #structure
        for k=1:length(T) #wavelets
            uwdot::Float64=0
        
            uwddot::Float64=(wav[k,1]-c[i].*uwdot-ks[i].*uw[1])/m
            
            for j=1:length(t)-1
                ph::Float64=wav[k,j+1]+a1[i]*uw[j]+a2[i]*uwdot+a3[i]*uwddot
                uw[j+1]=ph/kh[i]
                uwdoti::Float64=uwdot
                uwdot=gamma/(betha*increment).*(uw[j+1]-uw[j])+(1-gamma/betha).*uwdot+increment*(1-gamma/2/betha)*uwddot
                uwddot=1/(betha*increment^2).*(uw[j+1]-uw[j])-1/(betha*increment).*uwdoti-(1/2/betha-1)*uwddot
            end
            
            # Find peak time of each wavelet
            @. absltw::Vector{Float64}=abs(uw[:])*w[i]^2
            Rw::Float64=maximum(absltw)
            tw_index[k]=convert(Int64,findall(x->x==Rw,absltw)[1])
            tw_peak[k,i]=(tw_index[k]-1).*increment; 
            C[i,k]=uw[convert(Int64,t_index[i])].*w[i]^2;
        end 
    end
    # tw_peak=diag(tw_peak)
    @. Cini::Vector{Float64}=Cgain*C
    # Reduced off-diagonal C matrix
    DF::Vector{Int64} = setdiff(1:length(C),1:size(C)[1]+1:length(C))
    C[DF]=CoffDiag.*C[DF]; #suggested value 0.7
    # b matrix
    #println(C[end,:].*1000000)
    b=C\misfit

    atot::Vector{Float64}=amod
    for i=1:length(T)
        atot=atot+b[i].*wav[i,:]'
        wavtot=wavtot+b[i].*wav[i,:]
    end  

    for i=1:length(T) #structure
        fill!(0,uw)
        uwdot::Float64=0
        for j=1:length(t)-1
            ph::Float64=atot[j+1]+a1[i]*uw[j]+a2[i]*uwdot+a3[i]*uwddot
            uw[j+1]=ph/kh[i]
            uwdoti::Float64=uwdot
            uwdot=gamma/(betha*increment).*(uw[j+1]-uw[j])+(1-gamma/betha).*uwdot+increment*(1-gamma/2/betha)*uwddot
            uwddot=1/(betha*increment^2).*(uw[j+1]-uw[j])-1/(betha*increment).*uwdoti-(1/2/betha-1)*uwddot
        end
        
        # Find peak time of each wavelet   
        @. abslt::Vector{Float64}=abs(uw)*w[i]^2
        R::Float64=maximum(abslt)
        tw_index1::Int64= findall(x->x==R,abslt)[1]
        apeaknew::Float64=uw[tw_index1]*w[i]^2
    
        if apeaknew<0
            target[i]=-abs(target[i])
        end   
        misfitini[i]=target[i]-apeaknew; 
    end
    return b,Cini,misfitini,wavtot,t_peak
end