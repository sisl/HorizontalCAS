# Helper functions
sameSense(pra::Int, ra::Int) = mod(pra,2)==mod(ra,2)
strengthen(pra::Int, ra::Int) = ra>pra
weaken(pra::Int,ra::Int) = ra < pra
strongAlert(ra::Int) = ra > 2

# Reward function for HCAS MDP
function POMDPs.reward(mdp::HCAS_MDP, s::stateType, ra::actType)
    r = s[1]; t = s[2]; p = s[3]; vown = s[4]; vint = s[5]; pra = s[6]
    tau = mdp.currentTau
    rew = 0.0
    
    relx = r*cos(t); rely = r*sin(t); dx = vint*cos(p)-vown; dy = vint*sin(p)
    dv2 = dx*dx+dy*dy
    dCPA = r
    tCPA = 0.0
    if dv2>0.0
        tCPA = (-relx*dx - rely*dy)/dv2
        xT = relx + dx*tCPA
        yT = rely + dy*tCPA
        if tCPA>0.0
            dCPA = sqrt(xT*xT + yT*yT)
        else
            tCPA=0.0
        end
    end
    
    if tau==0.0
        if r<=500.0
            rew-=1.0
        else
            rew -= 1.0*exp(-(r-500.0)/500.0) # /300.0 for v5
        end
    end
    factor=1.0
    if pra!=COC
        factor=0.1 # factor always equals 1.0 for v5
    end
    if ra != COC
        rew-=5e-4 *factor#1e-9
        if strongAlert(ra)
            rew-=2e-3*factor
        end
        if (pra!=COC) && (!sameSense(pra,ra))
            rew-=5e-2
        elseif strengthen(pra,ra)
            rew-=1e-3
        elseif weaken(pra,ra)
            rew-=2e-4
        end
        
    # ra = COC
    else
        rew -= 1e-2*exp(-dCPA/500.0)*exp(-tCPA/10.0) /factor #dCPA/300.0 for v5
    end
    return rew
end