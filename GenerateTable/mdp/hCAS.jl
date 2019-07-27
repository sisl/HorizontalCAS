export HCAS_MDP

# Define MDP
mutable struct HCAS_MDP <: MDP{stateType, actType}
    ranges::Array{Float64,1}
    thetas::Array{Float64,1}
    psis::Array{Float64,1}
    vowns::Array{Float64,1}
    vints::Array{Float64,1}
    pras::Array{Int64,1}
    discount_factor::Float64
    turns::Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}
    currentTau::Float64
end
HCAS_MDP() = HCAS_MDP(RANGES,THETAS,PSIS,OWNSPEEDS,INTRSPEEDS,ACTIONS,discount_f,turns,0.0)
#HCAS_MDP(tau::Float64) = HCAS_MDP(RANGES,THETAS,PSIS,OWNSPEES,INTRSPEEDS,discount_f,turns,tau)

# Define necessary functions for HCAS MDP
POMDPs.actionindex(::HCAS_MDP, a::actType) = a + 1
POMDPs.actions(mdp::HCAS_MDP)  = ACTIONS
POMDPs.discount(mdp::HCAS_MDP) = mdp.discount_factor
POMDPs.n_actions(::HCAS_MDP)   = length(ACTIONS)

function POMDPs.convert_s(::Type{V} where V <: AbstractVector{Float64}, s::stateType, mdp::HCAS_MDP)
    v = [s[1],s[2],s[3],s[4],s[5],convert(Float64,s[6])]
    return v
end

function POMDPs.convert_s(::Type{stateType}, v::AbstractVector{Float64}, mdp::HCAS_MDP)
    s = (v[1],v[2],v[3],v[4],v[5],convert(Int,v[6]))
    return s
end