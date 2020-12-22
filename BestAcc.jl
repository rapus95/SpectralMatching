include("IniResponse.jl")
function BestAcc(amod::Vector{Float64},abest::Vector{Float64},t,m::Float64,tetha::Float64,dt::Float64,Tall::Vector{Float64},
                 wall::Vector{Float64},targetall::Vector{Float64},bestmisfitall::Float64)
    u,t_peak,t_index,apeak,misfit,_= IniResponse(amod,t,m,tetha,dt,Tall,wall,targetall);
    if bestmisfitall>=mean(abs.(misfit))
        bestmisfitall=mean(abs.(misfit)); 
        abest=amod; 
    end
    return abest,apeak,misfit,bestmisfitall
end