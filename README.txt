======================================
= README
======================================

Contents
1. Description of Data Set
2. Dynamic Linear Model

======================================
1. Description of Data Set
======================================

The available set is 107 points [DLMdata.csv]. 
As the amount of CO_2 changes points, it allows us to fit a 1st order dynamic linear model.

======================================
2. Dynamic Linear Model
======================================

Consider a time series (Y_t)t≥1 , in which we assume there is an un-observable Markov
chain (θ_t), called the state process, and Y_t is a measurement of θ_t containing error. Our
dynamic linear model could be specified as follows:
[1] y_{t} = F_{t} * θ_{t} + ν_{t}, ν_{t} ∼ N (0, V_{t})
[2] θ_{t} = G_{t} * θ_{t−1} + w_{t} , w_{t} ∼ N (0, W_{t})
for t = 1, ..., n where n = 107. We consider a Normal prior distribution for θ_{0} such that,
θ_{0} ∼ N (m_{0} , C_{0}),
in which m_{0} and C_{0} are known and equal to one. We follow the first order linear trend
mode and define:
V = σ^2 ,
W = σ_β^2 ,
followed by the notations proposed by Petris et al. (2007).
The unknown parameters are θ[states] and ψ = (σ^2 , σ_β^2 )'.

