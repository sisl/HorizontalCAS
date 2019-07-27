import h5py
import numpy as np

###### OPTIONS #######
tableFile = './Qtables/HCAS_oneSpeed_v6.h5'
trainingDataFiles = './TrainingData/HCAS_rect_TrainingData_v6_pra%d_tau%02d.h5'
######################

# Define state space. Make sure this matches up with the constants used to generate the MDP table!
acts = [0,1,2,3,4]
ranges = np.array([0.0,25.0,50.0,75.0,100.0,150.0,200.0,300.0,400.0,500.0,510.0,750.0,1000.0,1500.0,2000.0,3000.0,4000.0,5000.0,7000.0,9000.0,11000.0,13000.0,15000.0,17000.0,19000.0,21000.0,25000.0,30000.0,35000.0,40000.0,48000.0,56000.0])
thetas = np.linspace(-np.pi,np.pi,41)
psis  = np.linspace(-np.pi,np.pi,41)
vowns = [200.0]
vints = [200.0] 
taus  = np.linspace(0,60,61)

# Get table cutpoints
#X = np.array([[r,t,p,vo,vi] for vi in vints for vo in vowns for p in psis for t in thetas for r in ranges])

X = np.array([[r*np.cos(t),r*np.sin(t),p,] for p in psis for t in thetas for r in ranges])
#X = np.array([[r*np.cos(t),r*np.sin(t),p,vo,vi] for vi in vints for vo in vowns for p in psis for t in thetas for r in ranges])
# Compute means, ranges, mins and maxes
means = np.mean(X,axis=0)
rnges = np.max(X,axis=0)-np.min(X,axis=0)

# If only one value, then range is 0. Just divide by 1 instead of range
rnges = np.where(rnges==0.0,1.0,rnges)
X =(X-means)/rnges

# Compile table values
f = h5py.File(tableFile,'r')
Q = np.array(f['q'])
f.close()
Q = Q.T
ns2 = len(ranges)*len(thetas)*len(psis)*len(vowns)*len(vints)*len(acts)
ns3 = len(ranges)*len(thetas)*len(psis)*len(vowns)*len(vints)

print(Q.shape)
print(ns2*len(taus))
meanQ = np.mean(Q)
rangeQ = np.max(Q) - np.min(Q)
Q = (Q-meanQ)/rangeQ

means = np.concatenate((means,[meanQ]))
rnges = np.concatenate((rnges,[rangeQ]))

#min_inputs = np.array([ranges[0], thetas[0], psis[0], vowns[0], vints[0]])
#max_inputs = np.array([ranges[-1],thetas[-1],psis[-1],vowns[-1],vints[-1]])
min_inputs = np.array([-ranges[-1],-ranges[-1],psis[0]])
max_inputs = np.array([ ranges[-1], ranges[-1],psis[-1]])

#Save the Training Data
for tau in [0,5,10,15,20,30,40,60]:#range(100):
    Qsub = Q[tau*ns2:(tau+1)*ns2]
    for pra in range(5):
        Qsubsub = Qsub[pra*ns3:(pra+1)*ns3]
        with h5py.File(trainingDataFiles%(pra,tau),'w') as H:
            H.create_dataset('X',data=X)
            H.create_dataset('y',data=Qsubsub)
            H.create_dataset('means',data=means)
            H.create_dataset('ranges',data=rnges)
            H.create_dataset('min_inputs',data=min_inputs)
            H.create_dataset('max_inputs',data=max_inputs)
