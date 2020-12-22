function BroydenFunc(Broydeniter::Int64,dbgain::Float64,dCgain::Float64,bestmisfit::Float64,T::Array{Float64},t,
                    m::Float64,w::Array{Float64},tetha::Float64,increment::Float64,target::Array{Float64},b::Array{Float64},
                    u::Array{Float64},t_peak::Array{Float64},C::Array{Float64},misfit::Array{Float64},tol::Float64,
                    normtol::Float64,amod::Array{Float64},iter::Int64,wavmag::Float64,wavtot::Array{Float64},CoffDiag::Float64)

    a::Array{Float64}=amod
    gamma::Float64=0.5
    betha::Float64=0.25
    amod::Array{Float64}=a+wavtot'

    b=b'
    P::Int64=-1

    for q=1:Broydeniter #Broyden loop 
        unew::Array{Float64}=zeros(length(T),length(t))
        upeak_index::Array{Float64}=zeros(length(T))
        error::Array{Float64}=zeros(length(T))
        unew_peak::Array{Float64}=zeros(length(T))
        unewpeak_time::Array{Float64}=zeros(length(T))
        newmisfit::Array{Float64}=zeros(length(T),1)
        deltamisfit::Array{Float64}=zeros(length(T),1)
        wavtot::Array{Float64}=zeros(length(t)); #resetting total adjustment wavelets

        previousmeanmisfit::Float64=mean(abs.(misfit))
        Coff::Array{Float64}=P*C
        DF::Array{Int64} = setdiff(1:length(Coff),1:size(Coff)[1]+1:length(Coff))
        Coff[DF]=CoffDiag.*Coff[DF]; #suggested value 0.7
        deltab::Array{Float64}=dbgain.*(Coff\misfit)
        deltab = reshape(deltab,(size(deltab)[1],1))
        ns::Array{Float64}=deltab'*deltab
        b=b.+deltab'
        for j=1:length(T)    #wavelet
            uw::Array{Float64}=zeros(length(t))
            uwdot::Float64=0

            wj::Float64=w[j]*sqrt(1 .-tetha^2);  
            tshift::Float64=t_peak[j];#min(t_peak[j],3.9223*(1/Tn[j])^.-0.845)
            tj::Array{Float64}=t.-tshift;   
            gf::Float64=1.178.*(1/T[j])^.-0.93
            
            wav::Array{Float64}=wavmag.*b[j].*cos.(wj.*tj).*exp.(-(tj./gf).^2).*m
            ks::Float64=m.*w[j]^2
            c::Float64=2*m*w[j]*tetha
            uwddot::Float64=(wav[1].-c.*uwdot.-ks.*uw[1])/m
            
            a1::Float64=1/(betha*increment^2).*m+gamma/(betha*increment).*c
            a2::Float64=1/(betha*increment).*m+(gamma/betha-1).*c
            a3::Float64=(1/2/betha-1).*m+increment*(gamma/2/betha-1).*c
            kh::Float64=ks+a1
            
            for k=1:length(t)-1
                ph::Float64=wav[k+1]+a1*uw[k]+a2*uwdot+a3*uwddot
                uw[k+1]=ph/kh
                uwdoti::Float64=uwdot
                uwdot=gamma/(betha*increment).*(uw[k+1].-uw[k])+(1 .-gamma/betha).*uwdot+increment*(1 .-gamma/2/betha)*uwddot
                uwddot=1/(betha*increment^2).*(uw[k+1].-uw[k])-1/(betha*increment).*uwdoti.-(1/2/betha-1)*uwddot
            end
            
            #Find peak time of wavelet response
            abslt::Array{Float64}=abs.(uw);                      #Absolute wavelet pseudo acceleration response [psa]
            R::Float64=maximum(abslt);                       #Max of absolute psa
            tw_index::Int64=findall(x->x==R,abslt)[1];            #Index of max absolute psa on time series
            tw_peak::Float64=(tw_index-1).*increment;    #true value of total max peak     
            dTj::Float64=tw_peak.-t_peak[j];              #time shifting
            
            #Shifted wavelet
            tj=t.-tshift.+dTj; 
            wav=wavmag.*b[j].*cos.(wj.*tj).*exp.(-(tj./gf).^2).*m
            wavtot=wavtot.+wav';                 #Total wavelet
        end #end of j [wavelet]

        for i=1:length(T)    #structure
            uw::Array{Float64}=zeros(length(t))
            uwdot::Float64=0

            ks::Float64=m.*w[i]^2
            c::Float64=2*m*w[i]*tetha
            uwddot::Float64=(wavtot[1].-c.*uwdot.-ks.*uw[1])/m

            a1::Float64=1/(betha*increment^2).*m+gamma/(betha*increment).*c
            a2::Float64=1/(betha*increment).*m+(gamma/betha-1).*c
            a3::Float64=(1/2/betha-1).*m+increment*(gamma/2/betha-1).*c
            kh::Float64=ks+a1
        
            for j=1:length(t)-1
                ph::Float64=wavtot[j+1]+a1*uw[j]+a2*uwdot+a3*uwddot
                uw[j+1]=ph/kh
                uwdoti::Float64=uwdot
                uwdot=gamma/(betha*increment).*(uw[j+1].-uw[j])+(1 .-gamma/betha).*uwdot+increment*(1 .-gamma/2/betha)*uwddot
                uwddot=1/(betha*increment^2).*(uw[j+1].-uw[j])-1/(betha*increment).*uwdoti.-(1/2/betha-1)*uwddot
            end
        
            #Total response
            unew[i,:]=u[i,:]+uw
            
            #Calculate new misfit & deltamisfit
            abslt::Array{Float64}=abs.(unew[i,:])
            R::Float64=maximum(abslt)
            upeak_index[i]=findall(x->x==R,abslt)[1];
            unew_peak[i]=unew[i,convert(Int64,upeak_index[i])]; 
            unewpeak_time[i]=(upeak_index[i]-1)*increment
        
            if unew_peak[i]>=0
                target[i]=abs(target[i])
            else 
                target[i]=.-abs(target[i])
            end
    
            newmisfit[i]=target[i].-unew_peak[i]*w[i]^2; #New misfit
            deltamisfit[i]=newmisfit[i].-misfit[i];      #delta misfit
            misfit[i]=newmisfit[i];                     #Save new misfit as initial misfit for next iteration
            error[i]=maximum(abs.(target[i]).-abs.(unew_peak.*w[i]^2))/abs.(target[i])
        end # i loop [structure]
        #Broyden Updating
        C::Array{Float64}=C+dCgain.*((deltamisfit.-C*deltab).*deltab')./ns[1]
        meanmisfit::Float64=mean(abs.(misfit))
        maxerror::Float64=maximum(error)

        #println(iter," ------------ ", meanmisfit)
        #Save the best modified acceleration for each subperiod solved
        if bestmisfit>=meanmisfit
            bestmisfit=meanmisfit
            amod=a+wavtot'; 
        end
        iter=iter+1

        #Termination condition
        if meanmisfit>previousmeanmisfit
            P=.-P
        end; # I added this to aid solution convergence. In case that the misfit is diverging; the Broyden will change it direction
        if meanmisfit>1e3
            break
        end
        if meanmisfit<=normtol || maxerror<tol
            break
        end
    end
    return amod,bestmisfit,iter
end #end of broyden
