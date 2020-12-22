using DelimitedFiles
include("MatchingFunction.jl")
acc = readdlm("ElCentro.txt", '\t', Float64, '\n');
#Period Subset & Target
Tall=[0.05 0.075 0.1 0.12 0.15 0.175 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.5 3.0];                   #Input the entire period to be macthed in ascending order
targetall = [0.65,0.775,0.90,1,1,1,1,1,1,1,1,1,1,1,0.857142857142857,0.75,0.666666666666667,0.60,0.50,0.428571428571429,0.375000000000000,0.333333333333333,0.30,0.24,0.2]

#General Input
dt=0.02;               #Ground motion time step [sec]
tetha=0.05;             #Damping level
avgtol=0.01;            #Tolerance on average misfit
errortol=0.075;         #Tolerance on maximum error; 0.1 equals to 10#
zeropad=0;              #Length of zero pad at the beginningg & end of the ground motion [sec]

#Period Subset & Target
T1range=[0.07 1.0];     #Input the minimum & maximum period in Period subset 1. Example of valid use=[0.05 1.5]
T2range=[0.07 3.0];     #Input the minimum •& maximum period in Period subset 2  Example of valid use=[0.05 3.0], use 0 if no period subset 2
T3range=[0.0];            #Input the minimum & maximum period in Period subset 3, Example of valid use=[0.05 6.0]', use 0 if no period subset 3

#Gain Factors
Cgain=1e2;              #Gain factor on initial C
wavmag=1e-7;            #Wavelet magnitude; g  
CoffDiag=0.7;           #Coefficient of the off-diagonal terms  
dbgain=1.0;               #Gain factor on vector b updating
dCgain=1.0;               #Gain factor on C updating
Broydeniter=10;         #Maximum number of Broyden iteration
Outeriter=6;            #Maximum number of Outer-loop iteration
time = vcat(0:dt:(length(acc)-1)*dt);

#Execution
#@code_warntype SpectralMatchingFunction(acc,dt,tetha,avgtol,errortol,zeropad,Cgain,wavmag,CoffDiag,dbgain,dCgain,Broydeniter,Outeriter,Tall,T1range,T2range,T3range,targetall)
@time amod = SpectralMatchingFunction(acc,dt,tetha,avgtol,errortol,zeropad,Cgain,wavmag,CoffDiag,dbgain,dCgain,Broydeniter,Outeriter,Tall,T1range,T2range,T3range,targetall)
println(maximum(amod))
