function alpha_bounds(umin,umax)
    #= Upper and lower bounds for function e_1(u) = (u-sin(u))/u
       given upper/lower bounds on u
    Inputs:
        umin (array of floats): Lower bound on u for each cell
        umax (array of floats): Upper bound on u for each cell
    Outputs:
        (array of floats): Lower and upper bounds of error function
    =#
    num = length(umin)
    e1_umin = zeros(num); e1_umax = zeros(num)
    min_map = umin.!=0.0
    max_map = umax.!=0.0
    e1_umin[min_map] = 1 .-sin.(umin[min_map])./umin[min_map]
    e1_umax[max_map] = 1 .-sin.(umax[max_map])./umax[max_map]
    
    min_map = e1_umax .< e1_umin
    
    e1_min = copy(e1_umin)
    e1_max = copy(e1_umax)
    e1_max[min_map] = e1_umin[min_map]
    e1_min[min_map] = e1_umax[min_map]
    sign_map = (umin.<0.0) .& (umax.>0.0)
    e1_min[sign_map] .*=0.0
    return hcat(e1_min,e1_max)
end

function beta_bounds(umin,umax)
    #= Upper and lower bounds for function e_2(u) = (1-cos(u))/u
       given upper/lower bounds on u
    Inputs:
        umin (array of floats): Lower bound on u for each cell
        umax (array of floats): Upper bound on u for each cell
    Outputs:
        (array of floats): Lower and upper bounds of error function
    =#
    num = length(umin)
    e2_min = zeros(num); e2_max = zeros(num)
    min_map = umin.!=0.0
    max_map = umax.!=0.0
    
    e2_min[min_map] = (1 .-cos.(umin[min_map]))./umin[min_map]
    e2_max[max_map] = (1 .-cos.(umax[max_map]))./umax[max_map]
    return hcat(e2_min,e2_max)
end

function getNextSetBounds!(xBounds,yBounds,psiBounds,u0min,u0max,u1min,u1max;v0=200,v1=200)
    #= Get an over-approximated region of space reachable from each cell
    Inputs: 
        xBounds (array of floats): lower/upper bounds on x for each cell
        yBounds (array of floats): lower/upper bounds on y for each cell
        psiBounds (array of floats): lower/upper bounds on psi for each cell
        u0min (array of floats): lower bounds on u0 (ownship turn rate) for each cell
        u0max (array of floats): Upper bounds on u0 (ownship turn rate) for each cell
        u1min (array of floats): lower bounds on u1 (intruder turn rate) for each cell
        u1max (array of floats): Upper bounds on u1 (intruder turn rate) for each cell
        v0 (float): Ownship speed (ft/s)
        v1 (float): Intruder speed (ft/s)
    Outputs:
        None (xBounds, yBounds, and psiBounds are changed to be the boundaries of reachable space for each cell
    =#
    dpsi = Float32(0.5)*(psiBounds[:,2].-psiBounds[:,1])
    v1=Float32(v1)
    v0=Float32(v0)
    
    xBounds.+= v0.*alpha_bounds(u0min,u0max) .+ v1.*cos.((psiBounds[:,1].+psiBounds[:,2])./Float32(2.0)).-v0
    yBounds.-=v0*beta_bounds(u0min,u0max)[:,end:-1:1] .- v1.*sin.((psiBounds[:,1].+psiBounds[:,2])./Float32(2.0))
    
    v1cpsi = v1.*cos.((psiBounds[:,1].+psiBounds[:,2])./Float32(2.0))
    v1cpsi_pos = v1cpsi.>0.0
    v1cpsi_neg = .!v1cpsi_pos
    xBounds[v1cpsi_pos,:].-= v1cpsi[v1cpsi_pos].*alpha_bounds(u1min[v1cpsi_pos],u1max[v1cpsi_pos])[:,end:-1:1]
    xBounds[v1cpsi_neg,:].-= v1cpsi[v1cpsi_neg].*alpha_bounds(u1min[v1cpsi_neg],u1max[v1cpsi_neg])
    yBounds[v1cpsi_pos,:].+= v1cpsi[v1cpsi_pos].*beta_bounds(u1min[v1cpsi_pos],u1max[v1cpsi_pos])
    yBounds[v1cpsi_neg,:].+= v1cpsi[v1cpsi_neg].*beta_bounds(u1min[v1cpsi_neg],u1max[v1cpsi_neg])[:,end:-1:1]
    
    xBounds[v1cpsi_pos,1].+= v1cpsi[v1cpsi_pos].*(cos.(dpsi[v1cpsi_pos]).-1)
    yBounds[v1cpsi_pos,2].+= v1cpsi[v1cpsi_pos].*sin.(dpsi[v1cpsi_pos])
    yBounds[v1cpsi_pos,1].-= v1cpsi[v1cpsi_pos].*sin.(dpsi[v1cpsi_pos])    
        
    xBounds[v1cpsi_neg,2].+= v1cpsi[v1cpsi_neg].*(cos.(dpsi[v1cpsi_neg]).-1)
    yBounds[v1cpsi_neg,2].-= v1cpsi[v1cpsi_neg].*sin.(dpsi[v1cpsi_neg])
    yBounds[v1cpsi_neg,1].+= v1cpsi[v1cpsi_neg].*sin.(dpsi[v1cpsi_neg])
    v1cpsi = nothing; v1cpsi_pos=nothing; v1cpsi_neg=nothing; GC.gc()
    
    v1spsi = v1.*sin.((psiBounds[:,1].+psiBounds[:,2])./Float32(2.0))
    v1spsi_pos = v1spsi.>0.0
    v1spsi_neg = .!v1spsi_pos
    xBounds[v1spsi_pos,2].+= v1spsi[v1spsi_pos].*sin.(dpsi[v1spsi_pos])
    xBounds[v1spsi_pos,1].-= v1spsi[v1spsi_pos].*sin.(dpsi[v1spsi_pos])
    xBounds[v1spsi_pos,:].-= v1spsi[v1spsi_pos].*beta_bounds(u1min[v1spsi_pos],u1max[v1spsi_pos])[:,end:-1:1]
    yBounds[v1spsi_pos,1].+= v1spsi[v1spsi_pos].*(cos.(dpsi[v1spsi_pos]).-1)
    yBounds[v1spsi_pos,:].-= v1spsi[v1spsi_pos].*alpha_bounds(u1min[v1spsi_pos],u1max[v1spsi_pos])[:,end:-1:1]
    
    
    xBounds[v1spsi_neg,2].-= v1spsi[v1spsi_neg].*sin.(dpsi[v1spsi_neg])
    xBounds[v1spsi_neg,1].+= v1spsi[v1spsi_neg].*sin.(dpsi[v1spsi_neg])
    xBounds[v1spsi_neg,:].-= v1spsi[v1spsi_neg].*beta_bounds(u1min[v1spsi_neg],u1max[v1spsi_neg])
    yBounds[v1spsi_neg,2].+= v1spsi[v1spsi_neg].*(cos.(dpsi[v1spsi_neg]).-1)
    yBounds[v1spsi_neg,:].-= v1spsi[v1spsi_neg].*alpha_bounds(u1min[v1spsi_neg],u1max[v1spsi_neg])    
    v1spsi=nothing; v1spsi_pos=nothing; v1spsi_cos=nothing; GC.gc()
    
    psiBounds[:,2].+=u1max.-u0min
    psiBounds[:,1].+=u1min.-u0max
    
    cos_u0_min = cos.(max.(abs.(u0max),abs.(u0min)))
    sin_u0_min = sin.(u0min)
    sin_u0_max = sin.(u0max)
    
    xBounds_rotMin = min.(xBounds[:,1],cos_u0_min.*xBounds[:,1])+min.(sin_u0_min.*yBounds[:,1],sin_u0_max.*yBounds[:,1],sin_u0_min.*yBounds[:,2],sin_u0_max.*yBounds[:,2])
    xBounds_rotMax = max.(xBounds[:,2],cos_u0_min.*xBounds[:,2])+max.(sin_u0_min.*yBounds[:,1],sin_u0_max.*yBounds[:,1],sin_u0_min.*yBounds[:,2],sin_u0_max.*yBounds[:,2])
    
    yBounds_rotMin = min.(yBounds[:,1],cos_u0_min.*yBounds[:,1])-max.(sin_u0_min.*xBounds[:,1],sin_u0_max.*xBounds[:,1],sin_u0_min.*xBounds[:,2],sin_u0_max.*xBounds[:,2])
    yBounds_rotMax = max.(yBounds[:,2],cos_u0_min.*yBounds[:,2])-min.(sin_u0_min.*xBounds[:,1],sin_u0_max.*xBounds[:,1],sin_u0_min.*xBounds[:,2],sin_u0_max.*xBounds[:,2])
    
    xBounds[:,:] = hcat(xBounds_rotMin,xBounds_rotMax)
    yBounds[:,:] = hcat(yBounds_rotMin,yBounds_rotMax)
    cos_u0_min=nothing; sin_u0_min=nothing; sin_u0_max=nothing; xBounds_rotMin = nothing; xBounds_rotMax=nothing; yBounds_rotMin=nothing; yBounds_rotMax=nothing; GC.gc()
    
end

function getCutpoints(XS)
    xStepsAll = XS[2:end]-XS[1:end-1]
    xSteps = Array{Float32,1}()
    xNums = Array{Int32,1}()
    xCuts = [XS[1]]
    currentStep = -1
    for (i,step) in enumerate(xStepsAll)
        if step != currentStep
            currentStep = step
            xSteps = vcat(xSteps,[step])
            xNums = vcat(xNums,[i])
            if length(xNums)>1
                xCuts = vcat(xCuts,[xCuts[end]+(xNums[end]-xNums[end-1])*xSteps[end-1]])
            end
        end
    end
    return (xSteps, xCuts, xNums)
end

function getIndexBounds(xBoundsNext,yBoundsNext,psiBoundsNext)
    #= Compute which cells overlap with region of state space for each region
    Inputs:
        xBoundsNext (array of floats): Lower/upper bounds on x for each region
        yBoundsNext (array of floats): Lower/upper bounds on y for each region
        psiBoundsNext (array of floats): Lower/upper bounds on psi for each region
    Outputs:
        xInds (array of ints): First/last index of x-array that overlap with the region
        yInds (array of ints): First/last index of y-array that overlap with the region
        psiInds (array of ints): First/last index of psi-array that overlap with the region
    =#
    
    # Only use cell if it is inside, not just touching the boundary. This small adjustment 
    # prevents cells from being added that just touch the boundary of the region
    psiBoundsNext[:,1].+=1e-6
    psiBoundsNext[:,2].-=1e-6
    xBoundsNext[:,1].+=1e-6
    xBoundsNext[:,2].-=1e-6
    yBoundsNext[:,1].+=1e-6
    yBoundsNext[:,2].-=1e-6
        
    # The X/Y arrays of the grid are not uniform, but because it is piecewise linearly spaced, we
    # can compute the indices more efficiently than a brute force search
    xInds = ones(Int32,size(xBoundsNext)[1],2)
    xMap = xInds.>0
    xSteps, xCuts, xNums = getCutpoints(XS)
    for i=1:length(xSteps)
        xMap = xMap .& (xBoundsNext.>xCuts[i])
        xInds[xMap] = Int32.(div.(xBoundsNext[xMap].+(xSteps[i]*xNums[i]-xCuts[i]),xSteps[i]))
    end
    
    yInds = ones(Int32,size(yBoundsNext)[1],2)
    yMap = yInds.>0
    ySteps, yCuts, yNums = getCutpoints(YS)
    for i=1:length(ySteps)
        yMap = yMap .& (yBoundsNext.>yCuts[i])
        yInds[yMap] = Int32.(div.(yBoundsNext[yMap].+(ySteps[i]*yNums[i]-yCuts[i]),ySteps[i]))
    end
        
    # Psi is split uniformly, so mapping to indices is much quicker
    psiInds = ones(Int32,size(psiBoundsNext)[1],2)
    psiMap = psiInds.>0
    psiSteps, psiCuts, psiNums = getCutpoints(PSIS_DEG)
    for i=1:length(psiSteps)
        psiMap = psiMap .& (psiBoundsNext.>psiCuts[i])
        psiInds[psiMap] = Int32.(div.(psiBoundsNext[psiMap].*180.0./pi.+(psiSteps[i]*psiNums[i]-psiCuts[i]),psiSteps[i]))
    end
    
    xInds[xInds.>NUMX] .= NUMX
    yInds[yInds.>NUMY] .= NUMY
    psiInds[psiInds.>NUMP] .= NUMP
    
    psiInds[psiInds[:,1].>psiInds[:,2],1] .-=NUMP

    return xInds,yInds,psiInds
end


function getReachDynamics(xBoundInds,yBoundInds,psiBoundInds)
    #= Convert the first/last indices for each dimension to a sparse
       representation of the reachable dynamics. Output r has length the same as the number of cells.
       The cells reachable from cell i are given by c[r[i]:r[i+1]-1]
    
    This function is slow due to for-loops
    
    Inputs:
        xBoundInds (array of ints): First/last index of x-array that overlap with the region
        yBoundInds (array of ints): First/last index of y-array that overlap with the region
        psiBoundInds (array of ints): First/last index of psi-array that overlap with the region
    Outputs:
        r (array of ints): Pointers to the c array to determine where to begin reading
        c (array of ints): Array of cell indices
    =#
    dX = Int32(maximum(xBoundInds[:,2]-xBoundInds[:,1]))
    dY = Int32(maximum(yBoundInds[:,2]-yBoundInds[:,1]))
    dPsi = Int32(maximum(psiBoundInds[:,2]-psiBoundInds[:,1]))
    
    r = Vector{Int32}()
    c = Vector{Int32}()
    inds = Int32(1):NUMREGIONS
    
    mapX = BitArray(undef,NUMREGIONS)
    mapX .= true
    mapY = BitArray(undef,NUMREGIONS)
    mapPsi = BitArray(undef,NUMREGIONS)
    
    for i = Int32(0):dX
        println(@sprintf("%d of %d",i+1,dX+1))
        mapX .&= (xBoundInds[:,1].+i .<= xBoundInds[:,2])
        mapY.=mapX
        @time for j = Int32(0):dY
            mapY .&= (yBoundInds[:,1].+j .<= yBoundInds[:,2])
            mapPsi.=mapY
            for k = Int32(0):dPsi
                mapPsi .&= (psiBoundInds[:,1].+k .<= psiBoundInds[:,2])
                append!(c,inds[mapPsi])
                append!(r,indicesToIndex_Vector(psiBoundInds[mapPsi,1].+k,yBoundInds[mapPsi,1].+j,xBoundInds[mapPsi,1].+i))
            end
        end
    end
    sortIdx = sortperm(c)
    return r[sortIdx],c[sortIdx]
end