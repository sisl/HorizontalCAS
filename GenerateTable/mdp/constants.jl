export COC,WL,WR,SL,SR, stateType, actType, ACTIONS, discount_f, RANGES,THETAS,PSIS,OWNSPEEDS,INTRPSEEDS, interp, turns

# ADVISORY INDICES
COC=0
WL=1
WR=2
SL=3
SR=4

# State Type:
stateType = Tuple{Float64,Float64,Float64,Float64,Float64,Int}
actType = Int
ACTIONS = [COC,WL,WR,SL,SR]

# Default parameters
discount_f = 1.0

### STATE CUTPOINTS ###
RANGES = [0.0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,400.0,500.0,510.0,750.0,1000.0,1500.0,2000.0,3000.0,4000.0,5000.0,7000.0,9000.0,11000.0,13000.0,15000.0,17000.0,19000.0,21000.0,25000.0,30000.0,35000.0,40000.0,48000.0,56000.0]
THETAS = Array(LinRange(-pi,pi,41))
PSIS   = Array(LinRange(-pi,pi,41))
OWNSPEEDS = [200.0] #[100.0,200.0,300.0,400.0]  #ft/s
INTRSPEEDS = [200.0] #[100.0,200.0,300.0,400.0]  #ft/s

interp = LocalGIFunctionApproximator(RectangleGrid(RANGES,THETAS,PSIS,OWNSPEEDS,INTRSPEEDS,ACTIONS)) # Create the local function approximator using the grid


### Dictionaries to define transitions ###
probs = [0.5,0.25,0.25]
turns = Dict(COC=>([0.34,0.33,0.33],[0.0,1.5,-1.5].* pi/180.0 ),
              WL=>(probs,[1.5,2.0,1.25].* pi/180.0 ),
              WR=>(probs,[-1.5,-1.25,-2.0].* pi/180.0 ),
              SL=>(probs,[3.0,4.0,2.0].* pi/180.0 ),
              SR=>(probs,[-3.0,-2.0,-4.0].* pi/180.0 ),
              -1=>([0.34,0.33,0.33],[0.0,1.5,-1.5].* pi/180.0 )) # FOR v5, 0, 1, -1