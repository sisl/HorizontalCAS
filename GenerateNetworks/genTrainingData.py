import h5py
import numpy as np

###### OPTIONS #######
tableFile = './Qtables/HCAS_oneSpeed_v6.h5'
trainingDataFiles = './TrainingData/HCAS_rect_TrainingData_v6_pra%d_tau%02d.h5'
useRect = True      # Set to true to use rectangular coordinates (x,y) instead of polar coordinates(rho,theta) 
                    # as inputs to neural network. Changing to rectangular coordinates was used to ease reachability equations.
                    # The two coordinate systems are related using: x = rho * cos(theta), y = rho * sin(theta)
                    # 
oneSpeed = True     # If the speed dimension has only one value, set to True to remove the speed dimensions from network input
                    # If true, the neural network will have only three inputs: x, y, psi if useRect, or rho, theta, psi otherwise.
                    # If false, neural network will five inputs: x, y, psi, vown, vint if useRect, or rho, theta, psi, vown, vint otherwise
saveTaus = [0,5,10,15,20,30,40,60] # Tau values at which neural networks will be trained
######################

###############################
# Define state space. Make sure this matches up with the constants used to generate the MDP table!
acts = [0,1,2,3,4]
ranges = np.array([0.0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,400.0,500.0,510.0,750.0,1000.0,1500.0,2000.0,3000.0,4000.0,5000.0,7000.0,9000.0,11000.0,13000.0,15000.0,17000.0,19000.0,21000.0,25000.0,30000.0,35000.0,40000.0,48000.0,56000.0])
thetas = np.linspace(-np.pi,np.pi,41)
psis  = np.linspace(-np.pi,np.pi,41)
vowns = [200.0]
vints = [200.0] 
taus  = np.linspace(0,60,61)
###############################

# Get table cutpoints depending on useRect and oneSpeed settings
if useRect:
    if oneSpeed:
        X = np.array([[r,t,p] for p in psis for t in thetas for r in ranges])
    else:
        X = np.array([[r,t,p,vo,vi] for vi in vints for vo in vowns for p in psis for t in thetas for r in ranges])
else:
    if oneSpeed:
        X = np.array([[r*np.cos(t),r*np.sin(t),p,] for p in psis for t in thetas for r in ranges])
    else:
        X = np.array([[r*np.cos(t),r*np.sin(t),p,vo,vi] for vi in vints for vo in vowns for p in psis for t in thetas for r in ranges])

# Compute means, ranges, mins and maxes of inputs
means = np.mean(X, axis=0)
rnges = np.max(X, axis=0) - np.min(X, axis=0)
min_inputs = np.min(X, axis=0)
max_inputs = np.max(X, axis=0)

# Normalize each dimension of inputs to have 0 mean, unit range
# If only one value, then range is 0. Just divide by 1 instead of range
rnges = np.where(rnges==0.0, 1.0, rnges)
X  = (X - means) / rnges

# Load Q table
f = h5py.File(tableFile, 'r')
Q = np.array(f['q'])
f.close()
Q = Q.T

# Normalize entire output data to have 0 mean, unit range
# Add output normalization to the inputs means and rnges vector
meanQ = np.mean(Q)
rangeQ = np.max(Q) - np.min(Q)
Q = (Q - meanQ) / rangeQ
means = np.concatenate((means,[meanQ]))
rnges = np.concatenate((rnges,[rangeQ]))

# Sizes to help slice the table to create subtables used for training separate networks
ns2 = len(ranges) * len(thetas) * len(psis) * len(vowns) * len(vints) * len(acts)
ns3 = len(ranges) * len(thetas) * len(psis) * len(vowns) * len(vints)

#Save the Training Data
for tau in saveTaus:
    Qsub = Q[tau*ns2:(tau+1)*ns2]
    for pra in acts:
        Qsubsub = Qsub[pra*ns3:(pra+1)*ns3]
        with h5py.File(trainingDataFiles%(pra,tau),'w') as H:
            H.create_dataset('X', data=X)
            H.create_dataset('y', data=Qsubsub)
            H.create_dataset('means', data=means)
            H.create_dataset('ranges', data=rnges)
            H.create_dataset('min_inputs', data=min_inputs)
            H.create_dataset('max_inputs', data=max_inputs)
