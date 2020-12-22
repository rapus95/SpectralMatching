include("AccScale.jl")
include("BestAcc.jl")
include("BroydenFunc.jl")
include("IniResponse.jl")
include("InitialC2.jl")
include("PeriodSubset.jl")
include("Zeropad.jl")
using Statistics
function SpectralMatchingFunction(arec::Vector{Float64},dt::Float64,tetha::Float64,avgtol::Float64,errortol::Float64,zeropad::Int64,
                                  Cgain::Float64,wavmag::Float64,CoffDiag::Float64,dbgain::Float64,dCgain::Float64,Broydeniter::Int64,Outeriter::Int64,
                                  Tall::Vector{Float64},T1range::Vector{Float64},T2range::Vector{Float64},T3range,targetall::Vector{Float64})
    
    #Properties
    m::Float64=1;                        #mass
    iter::Int64=1;                     #starting iteration count
    wall::Vector{Float64}=2*pi./Tall

    #Zero padding at start & end of acceleration
    arec,t=Zeropad(zeropad,arec,dt)
    abest::Vector{Float64}=arec; 
    amod::Vector{Float64}=arec

    #Form the period subset
    T1::Vector{Float64},T2,T3,target1::Vector{Float64},target2,target3=PeriodSubset(Tall,T1range,T2range,T3range,targetall)
    
    #initial value of best misfit for all targets
    bestmisfitall::Float64=1000
    w::Vector{Float64}=2*pi./T1;
    @fastmath @inbounds for OuterLoop=1:Outeriter
        #Step 1; calculate initial response
        u,t_peak,t_index,apeak,misfit,target=IniResponse(amod,t,m,tetha,dt,T1,w,target1)
        #Step 2; Least square scalling
        amod=AccScale(apeak,target,amod)
        #Step 3: Scale the acceleration time series based on the least square fit
        u,t_peak,t_index,apeak,misfit,_=IniResponse(amod,t,m,tetha,dt,T1,w,target1)
        #Step 4; calculate initial b; initial C & initial misfit
        b,C,misfit,wavtot=InitialC2(T1,w,m,tetha,dt,t_peak,t_index,t,Cgain,target1,misfit,amod,CoffDiag,wavmag)
      
        #Step 5; Broyden Loop
        bestmisfit::Float64=mean(abs.(misfit))
        amod,bestmisfit,iter=BroydenFunc(Broydeniter,dbgain,dCgain,bestmisfit,T1,t,m,w,tetha,dt,target1,b,u,t_peak,C,misfit,errortol,avgtol,amod,iter,wavmag,wavtot,CoffDiag)
        #Step 6; Saving best result
        abest,apeak,misfit,bestmisfitall=BestAcc(amod,abest,t,m,tetha,dt,Tall,wall,targetall,bestmisfitall)
        
        #Step 7; Termination condition
        maxmisfit::Float64=abs(maximum((abs.(apeak)-abs.(targetall))./abs.(targetall)))

        if mean(abs.(misfit))<=avgtol||maxmisfit<errortol
            break
        end

        #PERIOD SUBSET 2
        if T2!=0
            w=2*pi./T2; #angular frequency

            #Step 1; calculate initial response
            u,t_peak,t_index,apeak,misfit,target=IniResponse(amod,t,m,tetha,dt,T2,w,target2)

            #Step 2; Least square scalling
            amod=AccScale(apeak,target,amod)

            #Step 3: Scale the acceleration time series based on the least square fit
            u,t_peak,t_index,apeak,misfit,_=IniResponse(amod,t,m,tetha,dt,T2,w,target2)

            #Step 4; calculate initial b; initial C & initial misfit
            b,C,misfit,wavtot=InitialC2(T2,w,m,tetha,dt,t_peak,t_index,t,Cgain,target2,misfit,amod,CoffDiag,wavmag)

            #Step 5; Broyden Loop
            bestmisfit=mean(abs.(misfit))
            amod,bestmisfit,iter=BroydenFunc(Broydeniter,dbgain,dCgain,bestmisfit,T2,t,m,w,tetha,
                dt,target2,b,u,t_peak,C,misfit,errortol,avgtol,amod,iter,wavmag,wavtot,CoffDiag)

            #Step 6; Saving best result
            abest,apeak,misfit,bestmisfitall=BestAcc(amod,abest,t,m,tetha,dt,Tall,wall,targetall,bestmisfitall)

            #Step 7; Termination condition
            maxmisfit=abs(maximum((abs.(apeak)-abs.(targetall))./abs.(targetall)))
            if mean(abs.(misfit))<=avgtol||maxmisfit<errortol
                break
            end
        end

        #PERIOD SUBSET 3
        if T3!=0; 
            w=2*pi./T3; #angular frequency

            #Step 1; calculate initial response
            u,t_peak,t_index,apeak,misfit,target=IniResponse(amod,t,m,tetha,dt,T3,w,target3)

            #Step 2; Least square scalling
            amod=AccScale(apeak,target,amod)

            #Step 3: Scale the acceleration time series based on the least square fit
            u,t_peak,t_index,apeak,misfit,_=IniResponse(amod,t,m,tetha,dt,T3,w,target3)

            #Step 4; calculate initial b; initial C & initial misfit
            b,C,misfit,wavtot=InitialC2(T3,w,m,tetha,dt,t_peak,t_index,t,Cgain,target3,misfit,amod,CoffDiag,wavmag)

            #Step 5; Broyden Loop
            bestmisfit=mean(abs.(misfit))
            amod,bestmisfit,iter=BroydenFunc(Broydeniter,dbgain,dCgain,bestmisfit,T3,t,m,w,tetha,
                dt,target3,b,u,t_peak,C,misfit,errortol,avgtol,amod,iter,wavmag,wavtot,CoffDiag)

            #Step 6; Saving best result
            abest,apeak,misfit,bestmisfitall=BestAcc(amod,abest,t,m,tetha,dt,Tall,wall,targetall,bestmisfitall)

            #Step 7; Termination condition
            maxmisfit=abs(maximum((abs.(apeak).-abs.(targetall))./abs.(targetall)))
            if mean(abs.(misfit))<=avgtol||maxmisfit<errortol
                break
            end
        end
    end
    return amod
end #of outer loop