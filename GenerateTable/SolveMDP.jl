@everywhere push!(LOAD_PATH,"mdp")
@everywhere using LocalFunctionApproximation
@everywhere using GridInterpolations
@everywhere using POMDPs
@everywhere using SparseArrays
@everywhere using HCAS
@everywhere using HDF5
@everywhere using POMDPModelTools
@everywhere using Printf


### OPTIONS ###
saveFile = "./Qtables/HCAS_oneSpeed_v6.h5"
nTau_warm=99
nTau_max = 60
###############

mdps = Array{Union{MDP, POMDP},1}()
for tau = 0:nTau_max
    push!(mdps, HCAS_MDP())
    mdps[end].currentTau = tau
end


@everywhere function compute_helper(states,n_states,mdps,grid,a,nTau_max)
    rval = zeros(Int32,n_states*100)
    cval = zeros(Int32,n_states*100)
    zval = zeros(Float32,n_states*100)
    rews = []
    for tau = 0:nTau_max
        push!(rews,zeros(n_states))
    end
    index=1
    fracBase=0.2
    frac = fracBase
    
    ## Iterate through all states
    for (i,s) in enumerate(states)
        if i/n_states>=frac
            print(round(frac*100))
            println("% Complete")
            frac+=fracBase
        end

        ## Compute rewards
        for tau=0:nTau_max
            rews[tau+1][i] = reward(mdps[tau+1],s,a)
        end

        ## Compute transitions
        dist = transition(mdps[1], s, a)
        for (sp, p) in weighted_iterator(dist)
            if p >0.0
                sp_point = POMDPs.convert_s(Vector{Float64}, sp, mdps[1])
                sps, probs = interpolants(grid, sp_point)
                for (spi, probi) in zip(sps,probs)
                    rval[index] = i
                    cval[index] = spi
                    zval[index] = probi*p
                    index += 1
                end
            end
        end
    end
    trans = sparse(rval[1:index-1],cval[1:index-1],zval[1:index-1],n_states,n_states)
    return (trans, rews)
end

@everywhere function compute_Qa(r,gam,trans,U)
    return r + gam*trans*U
end

function compute_trans_reward(mdps::Array{Union{MDP, POMDP},1},interp::LocalFunctionApproximator,nTau_max)
    ## Dictionaries to populate
    t = Dict()
    rews = Dict()
    rc = Dict()
    
    ## Compute states
    n_states = n_interpolating_points(interp)
    interp_points = get_all_interpolating_points(interp)
    S = statetype(typeof(mdps[1]))
    interp_states = Vector{S}(undef, n_states)
    for (i,pt) in enumerate(interp_points)
        interp_states[i] = POMDPs.convert_s(S, pt, mdps[1])
    end
    
    ## Loop through all actions in parallel
    for (ai,a) in enumerate(actions(mdps[1]))
        rc[ai] = remotecall(compute_helper, mod(ai,nprocs()-1)+2,interp_states,n_states,mdps,interp.grid,a,nTau_max)
    end
    for (ai,a) in enumerate(actions(mdps[1]))
        (t[ai], rews[ai]) = fetch(rc[ai])
    end
    return (t, rews)
end

function computeQ(mdps::Array{Union{MDP, POMDP},1},interp,nTau_warm,nTau_max)
    (trans, rews) = compute_trans_reward(mdps,interp,nTau_max);

    # Initialize Q and U vectors
    ns = length(rews[1][1])
    na = length(actions(mdps[1]))
    nt = nTau_max+1
    taus = 0:nTau_max
    gam = discount(mdps[1])
    U = zeros(ns)
    Q = zeros(ns,length(actions(mdps[1])))
    Q_out = zeros(ns*nt,na)
    Q_rc = Dict()

    #mdp_tau0 = []
    #mdp_tau1 = []

    ## Warm start for U at tau=0
    for i=1:nTau_warm
        @printf("Warm up: %d/%d\n",i,nTau_warm)
        for ai = 1:na
            Q_rc[ai] = remotecall(compute_Qa,mod(ai,nprocs()-1)+2,rews[ai][1],gam,trans[ai],U)
        end
        for ai = 1:na
            Q[:,ai] = fetch(Q_rc[ai])
        end
        U = maximum(Q,dims=2)
    end

    for i=1:nt
        @printf("Tau %d\n",taus[i])
        for ai = 1:na
            Q_rc[ai] = remotecall(compute_Qa,mod(ai,nprocs()-1)+2,rews[ai][i],gam,trans[ai],U)
        end
        for ai = 1:na
            Q[:,ai] = fetch(Q_rc[ai])
        end
        U = maximum(Q,dims=2)
        Q_out[1+(i-1)*ns:(i*ns),:] = deepcopy(Q)
    end
    return Q_out
end
@time Q_out = computeQ(mdps,interp,nTau_warm,nTau_max)


println("Writing Qvalues")
h5open(saveFile, "w") do file
    write(file, "q", Q_out) 
end
