function BroydenFunc(Broydeniter::Int64,dbgain::Float64,dCgain::Float64,bestmisfit::Float64,T::Vector{Float64},t,
                    m::Float64,w::Vector{Float64},tetha::Float64,increment::Float64,target::Vector{Float64},b::Vector{Float64},
                    u::Vector{Float64},t_peak::Vector{Float64},C::Vector{Float64},misfit::Vector{Float64},tol::Float64,
                    normtol::Float64,amod::Vector{Float64},iter::Int64,wavmag::Float64,wavtot::Vector{Float64},CoffDiag::Float64)

    a::Vector{Float64}=amod
    gamma::Float64=0.5
    betha::Float64=0.25
    amod::Vector{Float64}=a+wavtot'

    b=b'
    P::Int64=-1
    unew::Vector{Float64}=fill(0,(length(T),length(t)))
    upeak_index::Vector{Float64}=fill(0,length(T))
    error::Vector{Float64}=fill(0,length(T))
    unew_peak::Vector{Float64}=fill(0,length(T))
    unewpeak_time::Vector{Float64}=fill(0,length(T))
    newmisfit::Vector{Float64}=fill(0,length(T),1)
    deltamisfit::Vector{Float64}=fill(0,length(T),1)
    wavtot::Vector{Float64}=fill(0,length(t)); #resetting total adjustment wavelets
    uw::Vector{Float64}=Vector{Float64}(undef,length(t))
    @. tj::Vector{Float64} = t-tpeak
    @. Coff::Vector{Float64}=P*C
    DF::Vector{Int64} = setdiff(1:length(Coff),1:size(Coff)[1]+1:length(Coff))
    @. Coff[DF]=CoffDiag.*Coff[DF]; #suggested value 0.7
    @. ks::Vector{Int64}=m*w^2
    @. c::Vector{Int64}=2*m*w*tetha
    @. a1::Vector{Int64}=1/(betha*increment^2)*m+gamma/(betha*increment)*c
    @. a2::Vector{Int64}=1/(betha*increment)*m+(gamma/betha-1)*c
    @. a3::Vector{Int64}=(1/2/betha-1)*m+increment*(gamma/2/betha-1)*c
    @. kh::Vector{Int64}=ks+a1
    @. wj::Vector{Int64}=w*sqrt(1-tetha^2);  
    @. gf::Vector{Int64}=1.178*(1/T)^-0.93
    for q=1:Broydeniter #Broyden loop 
        previousmeanmisfit::Float64=mean(abs.(misfit))
        
        deltab::Vector{Float64}=dbgain.*(Coff\misfit)
        deltab = reshape(deltab,(size(deltab)[1],1))
        ns::Vector{Float64}=deltab'*deltab
        b=b.+deltab'
        
        fill!(0,uw)
        for j=1:length(T)    #wavelet          
            uwdot::Float64=0
            
            @. wav::Vector{Float64}=wavmag*b[j]*cos(wj[i]*tj[i])*exp((tj[i]/gf[i])^2)*m
            @. uwddot::Float64=(wav[1]-c[i]*uwdot-ks[i]*uw[1])/m
            
            for k=1:length(t)-1
                ph::Float64=wav[k+1]+a1[i]*uw[k]+a2[i]*uwdot+a3[i]*uwddot
                uw[k+1]=ph/kh[i]
                uwdoti::Float64=uwdot
                uwdot=gamma/(betha*increment).*(uw[k+1].-uw[k])+(1 .-gamma/betha).*uwdot+increment*(1 .-gamma/2/betha)*uwddot
                uwddot=1/(betha*increment^2).*(uw[k+1].-uw[k])-1/(betha*increment).*uwdoti.-(1/2/betha-1)*uwddot
            end
            
            #Find peak time of wavelet response
            abslt::Vector{Float64}=abs.(uw);                      #Absolute wavelet pseudo acceleration response [psa]
            R::Float64=maximum(abslt);                       #Max of absolute psa
            tw_index::Int64=findall(x->x==R,abslt)[1];            #Index of max absolute psa on time series
            tw_peak::Float64=(tw_index-1).*increment;    #true value of total max peak     
            dTj::Float64=tw_peak.-t_peak[j];              #time shifting
            
            #Shifted wavelet
            tj=t.-t_peak[j].+dTj; 
            wav=wavmag.*b[j].*cos.(wj[i].*tj[i]).*exp.(-(tj[i]./gf[i]).^2).*m
            wavtot=wavtot.+wav';                 #Total wavelet
        end #end of j [wavelet]

        for i=1:length(T)    #structure
            fill!(0,uw)
            uwdot::Float64=0

            uwddot::Float64=(wavtot[1].-c[i].*uwdot.-ks[i].*uw[1])/m

            for j=1:length(t)-1
                ph::Float64=wavtot[j+1]+a1[i]*uw[j]+a2[i]*uwdot+a3[i]*uwddot
                uw[j+1]=ph/kh[i]
                uwdoti::Float64=uwdot
                uwdot=gamma/(betha*increment).*(uw[j+1].-uw[j])+(1 .-gamma/betha).*uwdot+increment*(1 .-gamma/2/betha)*uwddot
                uwddot=1/(betha*increment^2).*(uw[j+1].-uw[j])-1/(betha*increment).*uwdoti.-(1/2/betha-1)*uwddot
            end
        
            #Total response
            unew[i,:]=u[i,:]+uw
            
            #Calculate new misfit & deltamisfit
            abslt::Vector{Float64}=abs.(unew[i,:])
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
        @. C::Vector{Float64}=C+dCgain*((deltamisfit-C*deltab)*deltab')/ns[1]
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
