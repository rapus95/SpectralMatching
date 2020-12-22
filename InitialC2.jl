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
    for i=1:length(T) #structure
        fill!(0,uw)
        uwdot::Float64=0

        wj::Float64=w[i]*sqrt(1-tetha^2);  
        tshift::Float64=t_peak[i]
        tj::Vector{Float64}=t.-tshift.+dTj[i];   
        gf::Float64=1.178.*(1/T[i])^-0.93
        wav[i,:]=@. wavmag*cos(wj*tj)*exp(-(tj/gf)^2)*m*b[i]
        #wavfunc!(wav,wavmag.*cos.(wj.*tj).*exp.(-(tj./gf).^2).*m.*b[i],i)
        
        ks::Float64=m.*w[i]^2
        c::Float64=2*m*w[i]*tetha
        uwddot::Float64=(wav[i,1]-c.*uwdot-ks.*uw[1])/m

        a1::Float64=1/(betha*increment^2).*m+gamma/(betha*increment).*c
        a2::Float64=1/(betha*increment).*m+(gamma/betha-1).*c
        a3::Float64=(1/2/betha-1).*m+increment*(gamma/2/betha-1).*c
        kh::Float64=ks+a1
        
        for j=1:length(t)-1
            ph::Float64=wav[i,j+1]+a1*uw[j]+a2*uwdot+a3*uwddot
            uw[j+1]=ph/kh
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
        tj=t.-tshift.+dTj[i]
        wav[i,:]=wavmag.*cos.(wj.*tj).*exp.(-(tj./gf).^2).*m.*b[i]
    end
   
    tw_peak=zeros(length(T),length(T))
    misfitini::Vector{Float64}=zeros(length(T))
    C::Vector{Float64}=ones(length(T),length(T))
    wavtot::Vector{Float64}=zeros(length(t))

    for i=1:length(T) #structure
        for k=1:length(T) #wavelets
            uwdot::Float64=0
        
            ks::Float64=m.*w[i]^2
            c::Float64=2*m*w[i]*tetha
            uwddot::Float64=(wav[k,1]-c.*uwdot-ks.*uw[1])/m

            a1::Float64=1/(betha*increment^2).*m+gamma/(betha*increment).*c
            a2::Float64=1/(betha*increment).*m+(gamma/betha-1).*c
            a3::Float64=(1/2/betha-1).*m+increment*(gamma/2/betha-1).*c
            kh::Float64=ks+a1
            
            for j=1:length(t)-1
                ph::Float64=wav[k,j+1]+a1*uw[j]+a2*uwdot+a3*uwddot
                uw[j+1]=ph/kh
                uwdoti::Float64=uwdot
                uwdot=gamma/(betha*increment).*(uw[j+1]-uw[j])+(1-gamma/betha).*uwdot+increment*(1-gamma/2/betha)*uwddot
                uwddot=1/(betha*increment^2).*(uw[j+1]-uw[j])-1/(betha*increment).*uwdoti-(1/2/betha-1)*uwddot
                #uddot[i,j+1]=1/(betha*increment^2).*(u[i,j+1]-u[i,j])-1/(betha*increment).*udot[i,j]-(1/2/betha-1)*uddot[i,j]
            end
            
            # Find peak time of each wavelet
            absltw::Vector{Float64}=abs.(uw[:]).*w[i]^2
            Rw::Float64=maximum(absltw)
            tw_index[k]=convert(Int64,findall(x->x==Rw,absltw)[1])
            tw_peak[k,i]=(tw_index[k]-1).*increment; 
            C[i,k]=uw[convert(Int64,t_index[i])].*w[i]^2;
        end 
    end
    # tw_peak=diag(tw_peak)
    Cini::Vector{Float64}=Cgain.*C
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
        ks::Float64=m.*w[i]^2
        c::Float64=2*m*w[i]*tetha
        uwddot::Float64=(atot[1]-c.*uwdot-ks.*uw[1])/m

        a1::Float64=1/(betha*increment^2).*m+gamma/(betha*increment).*c
        a2::Float64=1/(betha*increment).*m+(gamma/betha-1).*c
        a3::Float64=(1/2/betha-1).*m+increment*(gamma/2/betha-1).*c
        kh::Float64=ks+a1

        for j=1:length(t)-1
            ph::Float64=atot[j+1]+a1*uw[j]+a2*uwdot+a3*uwddot
            uw[j+1]=ph/kh
            uwdoti::Float64=uwdot
            uwdot=gamma/(betha*increment).*(uw[j+1]-uw[j])+(1-gamma/betha).*uwdot+increment*(1-gamma/2/betha)*uwddot
            uwddot=1/(betha*increment^2).*(uw[j+1]-uw[j])-1/(betha*increment).*uwdoti-(1/2/betha-1)*uwddot
        end
        
        # Find peak time of each wavelet   
        abslt::Vector{Float64}=abs.(uw).*w[i]^2
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