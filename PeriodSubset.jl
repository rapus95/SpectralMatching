function PeriodSubset(Tall::Vector{Float64},T1range::Vector{Float64},T2range::Vector{Float64},T3range::Vector{Float64},targetall::Vector{Float64})
    #This routine forms the period subset[s]

    #Period Subset1
    indexmin::Int64=argmin(abs.(Tall.-minimum(T1range)))[2];
    indexmax::Int64=argmin(abs.(Tall.-maximum(T1range)))[2];
    T1::Vector{Float64}=Tall[indexmin:indexmax];
    target1::Vector{Float64}=targetall[indexmin:indexmax];

    #Period Subset2
    if T2range[1]==0
        T2=0;
        target2==0
    else
        indexmin=argmin(abs.(Tall.-minimum(T2range)))[2];
        indexmax=argmin(abs.(Tall.-maximum(T2range)))[2];
        T2::Vector{Float64}=Tall[indexmin:indexmax];
        target2::Vector{Float64}=targetall[indexmin:indexmax];
    end

    #Period Subset3
    if T3range[1]==0
        T3=0
        target3=0
    else
        indexmin=argmin(abs.(Tall.-minimum(T3range)))[2];
        indexmax=argmin(abs.(Tall.-maximum(T3range)))[2];
        T3=Tall[indexmin:indexmax];
        target3=targetall[indexmin:indexmax];
    end
    return T1,T2,T3,target1,target2,target3
end