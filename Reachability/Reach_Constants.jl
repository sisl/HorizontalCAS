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
          LinRange(-1000,1000,51), LinRange(1040,6000,125),
          LinRange(6050,10000,80),  LinRange(10100,15000,50),
          LinRange(15200,30000,75), LinRange(31000,40000,10)))
YS = convert(Array{Float32,1},vcat(LinRange(-25000,-16000,10),LinRange(-15000,-8200,35),
          LinRange(-8000,-4100,40),LinRange(-4000,-1040,75),
          LinRange(-1000,1000,51),LinRange(1040,4000,75),
          LinRange(4100,8000,40),LinRange(8200,15000,35),LinRange(16000,25000,10)))

PSIS_DEG = convert(Array{Float32,1},LinRange(-180.0,180.0,181))
PSIS = convert(Array{Float32,1},PSIS_DEG*pi/180.0)

NUMX = Int32(length(XS)-1)
NUMY = Int32(length(YS)-1)
NUMP = Int32(length(PSIS)-1)
NUMACTION= Int32(length(ACTIONS))
NUMTAU = Int32(length(TAUS))
NUMREGIONS=NUMX*NUMY*NUMP
