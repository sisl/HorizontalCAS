# Advisory indices
COC=0
WL=1
WR=2
SL=3
SR=4

# Grid discretization
ACTIONS = [COC,WL,WR,SL,SR]
TAUS    = [0,5,10,15,20,30,40,60];

XS = convert(Array{Float32,1},vcat(LinRange(-10000,-5000,6),
          LinRange(-4000,-3250,4),LinRange(-3000,-1050,40),
          LinRange(-1000,1000,101), LinRange(1040,6000,125),
          LinRange(6050,10000,80),  LinRange(10100,15000,50),
          LinRange(15200,30000,75), LinRange(31000,40000,10)))
YS = convert(Array{Float32,1},vcat(LinRange(-25000,-16000,10),LinRange(-15000,-8200,35),
          LinRange(-8000,-4100,40),LinRange(-4000,-1040,75),
          LinRange(-1000,1000,101),LinRange(1040,4000,75),
          LinRange(4100,8000,40),LinRange(8200,15000,35),LinRange(16000,25000,10)))

PSIS = convert(Array{Float32,1},LinRange(-pi,pi,361))
NUMX = Int32(length(XS)-1)
NUMY = Int32(length(YS)-1)
NUMP = Int32(length(PSIS)-1)
NUMACTION= Int32(length(ACTIONS))
NUMTAU = Int32(length(TAUS))
NUMREGIONS=NUMX*NUMY*NUMP

# Numbers to help speed computation for grid discretizations
xSteps = [1000.0,250,50,20,40,50,100,200,1000]
xCuts = [-10000,-4000,-3000,-1000.0,1000,6000,10000,15000,30000]
xNums = [1, 7, 11, 51,151,276,356,406,481]

ySteps = [1000.0,200,100,40,20,40,100,200,1000]
yCuts = [-25000.0,-15000,-8000,-4000,-1000,1000,4000,8000,15000]
yNums = [1,11,46,86,161,261,336,376,411]