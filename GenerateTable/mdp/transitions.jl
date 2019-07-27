# State transition function
function POMDPs.transition(mdp::HCAS_MDP, s::stateType, ra::actType)
    r = s[1]; t = s[2]; p = s[3]; vown = s[4]; vint = s[5]; pra = s[6];
    
    # Computation is faster when using vector of static size
    nextStates = MVector{9, stateType}(undef)
    nextProbs = @MVector(zeros(9))
    next_pra = ra
    ind=1

    # Compute probabilities of next states using sigma point sampling
    ownProbs, ownTurns = mdp.turns[pra]
    intProbs, intTurns = mdp.turns[-1]
    for i = 1:3
        for j = 1:3
            next_r,next_t,next_p,next_vown,next_vint = dynamics(r,t,p,vown,vint,ownTurns[i],intTurns[j],pra)
            nextStates[ind] = (next_r,next_t,next_p,next_vown,next_vint,next_pra)
            nextProbs[ind]  = ownProbs[i]*intProbs[j]
            ind+=1
        end
    end

    return SparseCat(nextStates,nextProbs)
end

# Dynamic equations
function dynamics(r::Float64,t::Float64,p::Float64,vown::Float64,vint::Float64,ownTurn::Float64, intTurn::Float64, ra::Int)
    x_own = 0.0; y_own = 0.0; x_int = r*cos(t); y_int = r*sin(t)
    dx_own = vown; dy_own = 0.0; dx_int = vint*cos(p); dy_int = vint*sin(p)
    x_own += dx_own; y_own += dy_own; x_int += dx_int; y_int += dy_int
    
    x_int_new = x_int - x_own; y_int_new = y_int - y_own
    heading_own = ownTurn; heading_int = p + intTurn
    
    r_new = sqrt(x_int_new.^2 + y_int_new.^2)
    t_new = atan(y_int_new,x_int_new) - heading_own
    p_new = heading_int - heading_own
    
    t_new = t_new > pi  ? t_new-2*pi : t_new
    t_new = t_new < -pi ? t_new+2*pi : t_new
    p_new = p_new > pi  ? p_new-2*pi : p_new
    p_new = p_new < -pi ? p_new+2*pi : p_new
    
    return r_new, t_new, p_new, vown, vint
    
end