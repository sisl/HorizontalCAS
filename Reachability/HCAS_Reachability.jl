using Printf
using HDF5
using Mmap

include("NNet_Calculations.jl")
include("HCAS_Plotting.jl")
include("Reach_Functions.jl")
include("Reach_Constants.jl")


function getNetworks(folder="../networks",ver=2,hu=40,epochs=1000)
    #= Load all of the neural networks
    Given information about neural network version, load the networks into a list
    
    NOTE: File format assumed here that may need to be changed
    Inputs: 
        folder (string): Folder where neural networks are stored
        ver (int, 2): version number of networks
        hu  (int, 40): Number of hidden units in each layer (written in filename)
        epochs (int, 1000): Epoch number of neural network
    
    Outputs:
        nets: List of NNet objects
    =#
    nets = []
    for pra in ACTIONS
        nets_tau = []
        for tau in TAUS
            nnet_file_name = @sprintf("%s/HCAS_rect_v%d_pra%d_tau%02d_%dHU_%03d.nnet",folder,ver,pra,tau,hu,epochs)
            nets_tau = vcat(nets_tau,NNet(nnet_file_name))
        end
        if pra==COC
            nets = vcat(nets,nets_tau)
        else
            nets = hcat(nets,nets_tau)
        end
    end
    return nets    
end

function getAdvParams(ra,delta)
    #= Get a dictionary of parameters about the advisories
    u0 is the turn rate of the ownship, u1 is the turn rate of the intruder
    Inputs:
        ra (int): resolution advivisory (RA) index
        delta (float): value of delta to define acceleration bounds
    Outputs:
        params (dict): Dictionary of advisory parameters
    =#
    params=Dict()
    
    # Intruder turn rate bounds
    params["u1min"] = Float32(-delta *pi/180.0)
    params["u1max"] = Float32( delta *pi/180.0)
    
    # Ownship turn rate bounds (specific to RA)
    if ra==COC
        params["u0min"] = Float32( 0.0 *pi/180.0-delta *pi/180.0)
        params["u0max"] = Float32( 0.0 *pi/180.0+delta *pi/180.0)
    elseif ra==WL
        params["u0min"] =  Float32(1.5 *pi/180.0-delta *pi/180.0)
        params["u0max"] =  Float32(1.5 *pi/180.0+delta *pi/180.0)
    elseif ra==WR
        params["u0min"] =  Float32(-1.5 *pi/180.0-delta *pi/180.0)
        params["u0max"] =  Float32(-1.5 *pi/180.0+delta *pi/180.0)
    elseif ra==SL
        params["u0min"] =  Float32(3.5 *pi/180.0-delta *pi/180.0)
        params["u0max"] =  Float32(3.5 *pi/180.0+delta *pi/180.0)
    elseif ra==SR
        params["u0min"] =  Float32(-3.5 *pi/180.0-delta *pi/180.0)
        params["u0max"] =  Float32(-3.5 *pi/180.0+delta *pi/180.0)
    else
        print("Don't recognize this RA!")
    end
    return params
end

function getActionName(ra)
    #= Convert action index to a string name
    Inputs:
        ra (int): Advisory index
    Outputs:
        (string) name of advisory
    =#
    if ra == COC
        return "COC"
    elseif ra == WL
        return "WL"
    elseif ra == WR
        return "WR"
    elseif ra == SL
        return "SL"
    elseif ra == SR
        return "SR"
    else
        return "UNKNOWN RA"
    end
end

function indicesToIndex_Vector(pInds,yInds,xInds)
    #= Turn list of psi, y, and x indices into a vector of cell indexes
    Get the cell index for each set of psi, y, and x indices
    Inputs:
        pInds (list of int): List of psi index values
        yInds (list of int): List of y index values
        xInds (list of int): List of x index values
    Outputs:
        List of cell indexes
    =#
    yInds[yInds.<Int32(1)].+=NUMY
    pInds[pInds.<Int32(1)].+=NUMP
    return pInds .+ NUMP.*(yInds.-Int32(1)) .+ NUMP.*NUMY.*(xInds.-Int32(1))
end

function pointToIndex(x,y,psi)
    #= Determine the cell that includes the given point. Truncate values if needed
    Inputs:
        x (float): X-value of point
        y (float): Y-value of point
        psi (float): Psi-value of point
    Outputs:
        Cell index
    =#
    if psi>PSIS[end-1]
        psi=PSIS[end-1]
    end
    if y>YS[end-1]
        y=YS[end-1]
    end
    if x > XS[end-1]
        x=XS[end-1]
    end
    xInd = findall(XS.>x)[1]-1
    yInd = findall(YS.>y)[1]-1
    psiInd = findall(PSIS.>psi)[1]-1
    
    ind = indicesToIndex(psiInd,yInd,xInd)
    return ind
end

function pointToSet(x,y,psi)
    #= Given a point, return a one-hot vector of the cell indexes occupied by the point
    Inputs:
        x: X-value of point
        y: Y-value of point
        psi: Psi-value of point
    Outputs:
        One-hot vector of occupied cell
    =#
    ind = pointToIndex(x,y,psi)
    set = zeros(Int32,NUMREGIONS)
    set[ind] = Int32(1)
    return set
end

function indicesToIndex(pInd,yInd,xInd)
    #= Turn a single psi, y, and x index into a cell index
    Inputs:
        pInd (int): Psi index value
        yInd (int): Y index value
        xInd (int): X index value
    Outputs:
        Cell index
    =#
    if pInd <= NUMP && yInd <= NUMY && xInd <= NUMX && pInd>=1 && yInd>=1 && xInd>=1
        return pInd + NUMP*(yInd-1) + NUMP*NUMY*(xInd-1)
    end
    return -1
end

function indexToIndices(ind)
    #= Turn a cell index into psi, y, and x indices
    Inputs:
        Cell index (int)
    Outputs:
        pInd (int): Psi index value
        yInd (int): Y index value
        xInd (int): X index value
    =#
    pInd = mod(ind-1,NUMP)+1
    ind = div(ind-pInd,NUMP)
    yInd= mod(ind,NUMY)+1
    xInd = div(ind-yInd+1,NUMY)+1
    return (pInd,yInd,xInd)
end

function indexToBounds(ind)
    #= Turn a cell index into the hyper-rectangular bounds of the cell
    Inputs:
        Cell index (int)
    Outputs:
        psi Lower (float)
        psi Upper (float)
        y Lower (float)
        y Upper (float)
        x Lower (float)
        x Upper (float)
    =#
    pInd,yInd,xInd = indexToIndices(ind)
    return (PSIS[pInd],PSIS[pInd+1],YS[yInd],YS[yInd+1],XS[xInd],XS[xInd+1])
end

function getAllBounds()
    #= Turn a cell index into psi, y, and x indices
    Inputs:
        Cell index (int)
    Outputs:
        pInd (int): Psi index value
        yInd (int): Y index value
        xInd (int): X index value
    =#
    xMin = [XS[k]    for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    xMax = [XS[k+1]  for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    yMin = [YS[j]    for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    yMax = [YS[j+1 ] for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    pMin = [PSIS[i]     for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    pMax = [PSIS[i+1]   for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    return hcat(xMin,xMax), hcat(yMin,yMax), hcat(pMin,pMax)
end

function getTurnBounds(ra,delta)
    #= Create arrays for the turning min/max values for each region
    Each region has the same turn rate bounds in this case
    Inputs:
        act (int): action index
        delta: (float): parameter to relax turn rate bounds
    Outputs:
        u0min (float array): minimum turn rate of ownship for each region
        u0max (float array): maximum turn rate of ownship for each region
        u1min (float array): minimum turn rate of intruder for each region
        u1max (float array): maximum turn rate of intruder for each region
    =#
    params = getAdvParams(ra,delta)
    u0min = ones(Float32,NUMREGIONS)*params["u0min"]
    u0max = ones(Float32,NUMREGIONS)*params["u0max"]
    u1min = ones(Float32,NUMREGIONS)*params["u1min"]
    u1max = ones(Float32,NUMREGIONS)*params["u1max"]
    return u0min, u0max,u1min,u1max
end


######################################################
### Functions for different initial reachable sets ###
######################################################

function getInitialSet(;pd=0)
    #= Create an initial reachable set. The keys are the sequences of most recent 
       advisories. Initially, the only key is a sequence of COC, since the system
       has not taken effect yet for the initial reachable set. For this function,
       all cells are reachable initially.
    Inputs:
        pd (int): Seconds of pilot delay. This determines the number of most recent 
            advisories that we need to track.
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    reachSet = Dict()
    initialKey = string(COC)
    for i = 1:pd
        initialKey *= string(COC)
    end
    reachSet[initialKey] = BitArray(undef,NUMREGIONS)
    reachSet[initialKey].= true
    return reachSet
end

function getInitialSet_Outer(;pd=0)
    #= Same as getInitialSet(), except only cells far away horizontally are initialized 
       as reachable.
    Inputs:
        pd (int): Seconds of pilot delay. This determines the number of most recent 
            advisories that we need to track.
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    reachSet = Dict()
    initialKey = "0"
    for i = 1:pd
        initialKey *= "0"
    end
    reachSet[initialKey] = BitArray(undef,NUMREGIONS)
    reachSet[initialKey].= false
    regs = [((k==1) || (j==1) || (k==NUMX) || (j==NUMY)) for i=1:NUMP, j=1:NUMY, k=1:NUMX][:]
    reachSet[initialKey][regs].=true
    return reachSet
end

function getInitialSet_onTop(;pd=0,xbnd=4000.0,ybnd=4000.0)
    #= Same as getInitialSet(), except only cells with small horizontal distance are
       initialized as reachable.
    Inputs:
        pd (int): Seconds of pilot delay. This determines the number of most recent 
            advisories that we need to track.
        xbnd (float): Maximum distance in x-direction for a cell to be in reachable set
        ybnd (float): Maximum distance in y-direction for a cell to be in reachable set
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    reachSet = Dict()
    initialKey = "0"
    for i = 1:pd
        initialKey *= "0"
    end
    reachSet[initialKey] = BitArray(undef,NUMREGIONS)
    reachSet[initialKey].= false
    xBounds, yBounds, psiBounds = getAllBounds()
    regs = (xBounds[:,1].>-xbnd) .& (xBounds[:,2].<xbnd) .& (yBounds[:,1].>-ybnd) .& (yBounds[:,2].<ybnd)
    reachSet[initialKey][regs].=true
    return reachSet
end

function getInitialSet_point(pd=0;x=40000,y=0,psi=-pi)
    #= Same as getInitialSet(), except only the cell closest to the given point is
       included in the initial reachable set.
    Inputs:
        pd (int): Seconds of pilot delay. This determines the number of most recent 
            advisories that we need to track.
        x (float): X-value of point
        y (float): Y-value of point
        psi (float): Psi-value of point
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    reachSet = Dict()
    initialKey = "0"
    for i = 1:pd
        initialKey *= "0"
    end
    index = pointToIndex(x,y,psi)
    reachSet[initialKey] = BitArray(undef,NUMREGIONS)
    reachSet[initialKey].= false
    reachSet[initialKey][index] = true
    return reachSet
end


function getInitialSet_CPA(;time=10,minDist=2000.0,thickness=1000.0,v0 = 200.0)
    #= Same as getInitialSet(), except only cells that are in a ring around the NMAC
        region are kept. This ring is defined by the input parameters, and pilot delay
        is assumed to be 0 seconds.
    Inputs:
        v0 (float): Speed of aircraft (ft/s)
        time (float): Time of flight. Increasing makes the ring larger
        minDist (float): Minimum distance (ft) of intruder after aircraft travel at speed v0 for
            a time of t seconds. 
        thickness (float): Controls thickness of ring. Intruder must be at most 
            minDist+thickness (ft) from ownship to be included in reachable set
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    xBounds, yBounds, _ = getAllBounds()
    xMean = 0.5.*(xBounds[:,1].+xBounds[:,2]);
    yMean = 0.5.*(yBounds[:,1].+yBounds[:,2]);
    xBounds=nothing; yBounds=nothing; psiBounds=nothing; GC.gc()
    dist = sqrt.((xMean.-time*v0).^2 .+yMean.^2) .- time*v0
    regs = (dist.>=minDist) .& (dist.<=(minDist+thickness))
    
    keepReg = dist.<=(maxDist+thickness)
    nmacCells = sqrt.(xMean.^2 .+ yMean.^2).<520
    
    reachSet = Dict()
    for ra = ACTIONS
        initialKey = string(ra)
        reachSet[initialKey] = BitArray(undef,NUMREGIONS)
        reachSet[initialKey].= false
        reachSet[initialKey][regs].=true
    end
    return reachSet, keepReg, nmacCells
end

function getNmacCells()
    #= Compute cells that are NMACs (horizontal separation less than 500 feet)
    Inputs: 
        None
    Outputs:
        nmacCells (BitArray): True for cells that are in the NMAC region of the state space
    =#
    xBounds, yBounds, _ = getAllBounds()
    nmacCells = BitArray(undef,NUMREGIONS)
    nmacCells .= false
    
    # Focus search only on cells that could be within 500 foot radius
    potentialCells = (xBounds[:,2].>-500) .& (xBounds[:,1].<500) .& (yBounds[:,2].<500) .& (yBounds[:,1].<500)
    
    # Check each corner of cell to see if it falls within NMAC region (500 feet horizontally)
    for i=[1,2]
        for j=[1,2]
            nmacCells[potentialCells] .|= sqrt.(xBounds[potentialCells,i].^2 .+ yBounds[potentialCells,j].^2).<500
        end
    end
    return nmacCells
end
    

###############################################
### Helper functions for the reachable sets ###
###############################################

function saveSets(filename, sets)
    #= Write the sets dictionary to an h5py file
       The first key is the time
       The second key is the advisory sequence
    Inputs: 
        filename (string): Name of file to write
        sets (Dictionary): Reachable set dictionary
    =#
    h5open(filename,"w") do file
        for (k,v) in sets
            println(k)
            for (k2,v2) in v
                write(file,string(k)*"/"*k2,convert(Array{UInt8},v2))
            end
        end
    end
end

function pruneSets!(sets,regKeep)
    #= Removes cells from the reachable set that are not in regKeep
    Inputs: 
        sets (Dictionary): Reachable set. The key is the sequence of previous advisories
        regKeep (BitArray): Array that is true if we should keep that cell in the set
    Outputs:
        None (The sets input is modified)
    =#
    for (k,v) = sets
        v[.!regKeep].=false
    end
end

function copySet(set)
    #= Deep copy a dictionary of BitArrays
    Inputs:
        set (Dictionary of BitArrays): Reachable set
    Outputs:
        out: Copy of reachable set
    =#
    out = Dict()
    for (k,v) = set
        out[k] = copy(v)
    end
    return out
end

function isConverged(sets,t;verbose=false)
    #= Check if a set of cells has converged at time t. Convergence
        is checked by seeing if any change has been made at the previous step.
        Since t counts down with each step, check if the reachable set at t+1
        is identical to the reachable set at time t.
    Inputs:
        Sets (Dictionary): Reachable sets dictionary. The first key is tau, 
            second key is most recent advisory sequence
        t (int): Tau value to check for convergence
    Outputs:
        allTrue (bool): True if the reachable set has not changed since the previous step
    
    =#
    allTrue = true
    sumDiff = 0
    if !(t+1 in keys(sets))
        return false
    end
    for key in keys(sets[t])
        if key in keys(sets[t+1])
            allTrue &= sets[t][key]==sets[t+1][key]
            if verbose
                sumDiff += sum(sets[t][key].!=sets[t+1][key])
            elseif !allTrue
                return false
            end
        else
            allTrue = false
        end
    end
    if verbose
        @printf("Number of different cells: %d\n",sumDiff)
    end
    return allTrue
end

function isNmac(sets,t,nmacCells)
    #= Check if any cells are reachable that are NMACS
       (near midair collisions)
    Inputs:
        Sets (Dictionary): Reachable sets dictionary. The first key is tau, 
            second key is most recent advisory sequence
        t (int): Tau value to check
        nmacCells (BitArray): True for cells that are NMACs
    Outputs:
        (bool): True if any cells in set are NMAC cells
    =#
    for (k,v) = sets[t]
        if any(v .& nmacCells)
            return true
        end
    end
    return false
end


##############################################################################
### Functions to create a memory-mapped array for over-approximated policy ###
##############################################################################

function sampleToMmap(nets;vown=200.0,vint=200.0,batch_size = 10000000)
    #= Evaluate network by sampling points at the center of each cell
    Inputs:
        nets (List of NNet): Neural networks to evaluate
        vown (float): Ownship speed (ft/s)
        vint (float): Intruder speed (ft/s)
        batch_size (int): Number of points to evaluate at a time
    Outputs:
        acts (BitArray): True where the network could give an advisory in a cell
    =#
    acts = BitArray(undef,NUMTAU,NUMACTION,NUMACTION,NUMREGIONS)
    acts.=false
    
    # Check center of cell. Change this to check a different point within each cell
    dx = 0.5; dy = 0.5; dp = 0.5
    
    # Making Mesh
    println("Making Mesh")
    pCenter = (1-dp)*PSIS[2:end]+dp*PSIS[1:end-1]
    yCenter = (1-dy)*YS[2:end]+dy*YS[1:end-1]
    xCenter = (1-dx)*XS[2:end]+dx*XS[1:end-1]
    xMesh = [x for p=pCenter, y=yCenter, x=xCenter]
    yMesh = [y for p=pCenter, y=yCenter, x=xCenter]
    pMesh = [p for p=pCenter, y=yCenter, x=xCenter]
    vintMesh = ones(NUMREGIONS)*vint 
    vownMesh = ones(NUMREGIONS)*vown
    netIn = permutedims(hcat(xMesh[:],yMesh[:],pMesh[:],vownMesh[:],vintMesh[:]),[2,1])
    
    # Evaluate each network through batches
    # The output associated with the highest score is the advisory given by the network
    for (tauInd,tau) = enumerate(TAUS)
        for (praInd,pra) = enumerate(ACTIONS)
            @printf("Tau: %d, pRA: %d\n",tau,pra)
            net = nets[tauInd,praInd]
            currentInd = 1
            while currentInd<NUMREGIONS
                println(currentInd," of ",NUMREGIONS)
                nextInd = min(currentInd+batch_size,NUMREGIONS+1)
                netOut = evaluate_network_multiple(net,netIn[:,currentInd:nextInd-1])
                _, mxindx = findmax(netOut,dims=1)
                for i = 1:length(mxindx)
                    acts[tauInd,praInd,mxindx[i][1],currentInd+i-1] = true
                end
                currentInd=nextInd
            end
        end
    end
    return acts
end

function reluValToMmap(folder,ver=6;praInds=1:NUMACTION,tauInds=1:NUMTAU)
    #= Read the output from ReluVal to create the network action BitArray
    Inputs:
        ver (int): Version of neural networks
        praInds (Array of int): previous advisories to read
        tauInds (Array of int): tau indices to read
    Outputs:
        acts (BitArray): True where the network could give an advisory in a cell
    =#
    acts = BitArray(undef,NUMTAU,NUMACTION,NUMACTION,NUMREGIONS)
    acts .= false
    for praInd = praInds
        for tauInd = tauInds
            println(@sprintf("RA index: %d, Tau Index: %d\n",praInd,tauInd))
            file = @sprintf("%s/HCAS_v%d_pra%d_tau%02d.txt",folder,ver,ACTIONS[praInd],TAUS[tauInd])
            raInd=0
            f=open(file)
            while !eof(f)
                line=readline(f)
                if length(line)>0 && line[1]!='t' # Skip empty lines and time statements
                    if line[1]=='R' # If the first letter is an 'R', advance to next ra index
                        raInd+=1
                        println(raInd)
                    else
                        reg = parse(Int32,split(line,",")[1])+1 # +1 because Julia is 1-indexed
                        acts[tauInd,praInd,raInd,reg] = true; # Mark the tau/pra/ra/cell as true
                    end
                end
            end
            close(f)
        end
    end
    return acts
end


function writeNetworkActionsMmap(;folder="NetworkApprox",nnetFolder="../networks",reluvalFolder="../ReluVal/Results",ver=6, hu=25, epochs=3000,useReluVal=false)
    #= Write advisories given by network to a memory-mapped array. Advisories can
       be determined by sampling a point within each cell or by using ReluVal. If 
       using ReluVal, then ReluVal must be run first.
    Inputs:
        folder (string): Folder to write memory-mapped file
        ver (int): Version of neural network
        hu (int): Hidden units in each layer of network
        epochs (int): Number of epochs used for training
        useReluVal (bool): True if using ReluVal, False for sampling
    Outputs:
        None (results are written to memory-mapped files)
    =#
    
    token="samp"
    if useReluVal
        token="reluVal"
    end
    
    filename=@sprintf("%s/netActions_mm_v%d_%s_rect.bin",folder,ver,token)
    netActions=nothing
    if useReluVal
        @time netActions = reluValToMmap(reluvalFolder,ver)
    else
        nets = getNetworks(nnetFolder,ver,hu,epochs)
        @time netActions = sampleToMmap(nets);
    end
    
    # Write to memory mapped file. First write the dimensions of the array, then write the full array
    # Dimension of netActions are tau, previous advisory, next advisory, and cell
    # The cells are 3D and have dimensions x, y, and psi
    s = open(filename,"w+")
    write(s,Int64(size(netActions,1)))
    write(s,Int64(size(netActions,2)))
    write(s,Int64(size(netActions,3)))
    write(s,Int64(size(netActions,4)))
    write(s,netActions)
    close(s)
end

function readNetworkActionsMmap(;folder="NetworkApprox",ver=6,token="reluVal")
    #= Read advisories given by network from a memory-mapped array.
    Inputs:
        folder (string): Folder to write memory-mapped file
        ver (int): Version of neural network
        token (string): Indentifying string (samp or reluVal)
    Outputs:
        netActions: Memory-mapped 4-D BitArray
    =#
    filename=@sprintf("%s/netActions_mm_v%d_%s_rect.bin",folder,ver,token)
    s = open(filename)
    m1 = read(s,Int64)
    m2 = read(s,Int64)
    m3 = read(s,Int64)
    m4 = read(s,Int64)
    return Mmap.mmap(s,BitArray,(m1,m2,m3,m4))
end

###########################################################################
### Generate memory-mapped array of cell transitions as a sparse array  ###
###########################################################################

function writeReachDynamicsMmap(;folder="ReachDynamics",ras=ACTIONS,ver=1,delta=0.0,v0=200.0,v1=180.0)
    #= Write reachableDynamics for each advisory to a memory-mapped array
    Inputs:
        folder (string): Folder to write memory-mapped file
        ras (list of int): List of advisories to consider
        delta (float): Parameter to relax the turn rate bounds
        v0 (float): Speed of ownship (ft/s)
        v1 (float): Speed of intruder (ft/s)
    Outputs:
        None (results are written to memory-mapped files)
    =#
    for ra = ras
        @printf("RA: %s\n",getActionName(ra))
        r,c = computeReachDynamics(delta,ra,v0=v0,v1=v1)
        idx_ptr = [UInt32(1)]
        append!(idx_ptr,findall(c[1:end-1].!=c[2:end]).+UInt32(1))
        append!(idx_ptr,[UInt32(length(c)+1)])
        next_idx = UInt32.(r)

        s = open(@sprintf("%s/reachDynamics_mm_v%d_vown%.1f_vint%.1f_delta%.2f_ra%d_idxPtr.bin",folder,ver,v0,v1,delta,ra),"w+")
        write(s, UInt32(size(idx_ptr,1)))
        write(s, idx_ptr)
        close(s)

        s = open(@sprintf("%s/reachDynamics_mm_v%d_vown%.1f_vint%.1f_delta%.2f_ra%d_nextIdx.bin",folder,ver,v0,v1,delta,ra),"w+")
        write(s, UInt32(size(next_idx,1)))
        write(s, next_idx)
        close(s)
        r=nothing; c=nothing; idx_ptr=nothing; next_idx=nothing
    end
end

function readReachDynamicsMmap(;folder="ReachDynamics",ver=1,delta=0.0,v0=200.0,v1=180,ras=ACTIONS)
    #= Read reachableDynamics for each advisory as a memory-mapped array
    Inputs:
        folder (string): Folder to read memory-mapped file
        ver (int): Version of dynamics
        delta (float): Parameter to relax the turn rate bounds
        v0 (float): Speed of ownship (ft/s)
        v1 (float): Speed of intruder (ft/s)
        ras (list of int): List of advisories to load
    Outputs:
        reachDynamics: List of memory-mapped vectors
    =#
    reachDynamics = Array{Vector{UInt32},1}[]
    
    for ra = ras
        fn = @sprintf("%s/reachDynamics_mm_v%d_vown%.1f_vint%.1f_delta%.2f_ra%d_idxPtr.bin",folder,ver,v0,v1,delta,ra)
        s = open(fn)
        m = read(s,UInt32)
        idxPtr = Mmap.mmap(s,Vector{UInt32},m)

        fn = @sprintf("%s/reachDynamics_mm_v%d_vown%.1f_vint%.1f_delta%.2f_ra%d_nextIdx.bin",folder,ver,v0,v1,delta,ra)
        s2 = open(fn)
        m = read(s2,UInt32)
        nextIdx = Mmap.mmap(s2,Array{UInt32,1},(m,));
    
        push!(reachDynamics,[idxPtr,nextIdx])
    end
    return reachDynamics
    
end

##########################################################
### Functions for computing the reachable set dynamics ###
##########################################################

function computeReachDynamics(delta,ra;v0=200.0,v1=200.0)
    #= Compute the reachable set dynamics
    Inputs:
        delta (float): Paramter to relax turn rate bounds
        ra (int): Advisory
        v0 (float): Ownship speed (ft/s)
        v1 (float): Intruder speed (ft/s)
    Outputs:
        reachDynamics (list): Data structure representing which cells
            can be reached at the next time step from the current cell
    =#
    
    # Get the initial cell and turn rate bounds
    xBounds,yBounds,psiBounds = getAllBounds()
    u0min,u0max,u1min,u1max= getTurnBounds(ra,delta)
    
    # Compute the region of the state space reachable at the next time step
    getNextSetBounds!(xBounds,yBounds,psiBounds,u0min,u0max,u1min,u1max,v0=v0,v1=v1)
    
    # Clear some memory
    u0min=nothing;u0max=nothing;u1min=nothing;u1max=nothing; GC.gc()
    
    # Compute the cells that overlap with the reachable region
    xBoundInds, yBoundInds,psiBoundInds = getIndexBounds(xBounds,yBounds,psiBounds)
    
    # Clear some memory
    xBounds=nothing; yBounds=nothing; psiBounds=nothing; GC.gc()
    
    # Return reachable dynamics array represented with vectors of pointers and indices
    return getReachDynamics(xBoundInds,yBoundInds,psiBoundInds)
end

############################################################################
### Function to compute the set of reachable cells at the next time step ###
############################################################################

function getNextSetMmap(set,networkActions,reachDynamics,tau;alwaysTrueRegions=nothing)
    #= Compute the next reachable set
    Inputs:
        set (Dictionary of BitArrays): Current reachable set
        networkActions (BitArray): Memory-mapped array that is True where the network 
            could give the advisory in that cell.
        reachDynamics (List): Memory-mapped data structure used to determine which cells are
            reachable at the next time step give a current cell and acting advisory
        tau: Tau value, used to determine which network is being used 
        alwaysTrueRegions (BitArray): True for a cell that should always be reachable, such
            as cells that are far away horizontally. An intruder could always appear on the horizon
    Outputs:
        nextSet (Dictionary of BitArrays): Next reachable set
    =#
    
    # Determine the index of tau so we use the correct neural network
    if tau<TAUS[1]
        tau=TAUS[1]
    end
    tauInd = findall(tau.>=TAUS)[end]
    
    nextSet = Dict()
    kl=1 #length of key
    for (key,currentCells) in set
        kl = length(key)
        
        delayAdvs = [parse(Int32,key[i]) for i=1:length(key)]   # Advisory the pilot could be following, delayed
        pRA       = parse(Int32,key[end]) # Advisory just given by the network at the previous time step
    
        # Loop through each possible next advisory
        for ra=ACTIONS 
            
            # Get the next key and determine in which cells the advisory could be given
            cellsGiven = currentCells .& networkActions[tauInd,pRA+1,ra+1,:]
            if any(cellsGiven)
                nextKey = key[2:end]*string(ra)
                
                # Initialize reachable array if key is not yet in the dictionary
                if !(nextKey in keys(nextSet))
                    nextSet[nextKey] = BitArray(undef,NUMREGIONS)
                    nextSet[nextKey].=false
                end

                # Use reachDynamics to determine add cells to the next reachable set
                cg = findall(cellsGiven)
                for actingRA in delayAdvs
                    idxPtrs = reachDynamics[actingRA+1][1]
                    nextIdxs = reachDynamics[actingRA+1][2]
                    ns = nextSet[nextKey]

                    for cell = cg
                        ns[nextIdxs[idxPtrs[cell]:(idxPtrs[cell+1]-1)]].=true
                    end
                end
            end
        end
    end
    
    # Add cells that are always reachable
    if alwaysTrueRegions != nothing
        key = ""
        for i=1:kl
            key*=string(COC)
        end
        if !(key in keys(nextSet))
            nextSet[key] = BitArray(undef,NUMREGIONS)
            nextSet[key].= false
        end
        nextSet[key][alwaysTrueRegions] .= true   
    end
    return nextSet
end
