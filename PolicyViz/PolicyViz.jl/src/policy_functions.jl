export get_belief,get_qval!,Policy,read_policy,evaluate

mutable struct Policy
    alpha       :: Matrix{Float64}
    actions     :: Array{Int64,2}
    nactions    :: Int64
    qvals       :: Vector{Float64}
end

function Policy(alpha::Matrix{Float64}, actions::Array{Int64,2})
    return Policy(alpha, actions, size(actions,2), zeros(size(actions,2)))
end

function read_policy(actions::Array{Int64,2}, alpha::Matrix{Float64})
    return Policy(alpha, actions)
end # function read_policy

function evaluate(policy::Policy, belief::SparseMatrixCSC{Float64,Int64})
    fill!(policy.qvals, 0.0)
    get_qval!(policy, belief)
    return copy(policy.qvals)
end # function evaluate

function get_qval!(policy::Policy, belief::SparseMatrixCSC{Float64, Int64})
    fill!(policy.qvals, 0.0)
    for iaction in 1:policy.nactions
        for ib in 1:length(belief.rowval)
            policy.qvals[iaction] += belief.nzval[ib] * policy.alpha[belief.rowval[ib], iaction]
        end # for b
    end # for iaction
end # function get_qval!

function get_belief(pstate::Vector{Float64}, grid::RectangleGrid,interp::Bool=false)
    belief = spzeros(NSTATES, 1)
    indices, weights = interpolants(grid, pstate)
    if !interp
        indices = indices[findmax(weights)[2]]
        weights = 1.0
    end
    for i = 1:length(indices)
        belief[indices[i]] = weights[i]
    end # for i
    return belief
end # function get_belief