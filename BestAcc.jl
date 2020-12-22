include("IniResponse.jl")
function BestAcc(amod::Array{Float64},abest::Array{Float64},t,m::Float64,tetha::Float64,dt::Float64,Tall::Array{Float64},
                 wall::Array{Float64},targetall::Array{Float64},bestmisfitall::Float64)
    u,t_peak,t_index,apeak,misfit,_= IniResponse(amod,t,m,tetha,dt,Tall,wall,targetall);
    if bestmisfitall>=mean(abs.(misfit))
        bestmisfitall=mean(abs.(misfit)); 
        abest=amod; 
    end
    return abest,apeak,misfit,bestmisfitall
end