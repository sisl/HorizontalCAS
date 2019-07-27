import numpy as np
import math
import tensorflow as tf
import h5py
import sys
import os


######## OPTIONS #########
ver = 6              # Neural network version
table_ver = 6        # Table Version
hu = 25              # Number of hidden units in each hidden layer in network
saveEvery = 1000     # Epoch frequency of saving
totalEpochs = 3000   # Total number of training epochs
trainingDataFiles = "./TrainingData/HCAS_rect_TrainingData_v%d_pra%d_tau%02d.h5" # File format for training data
nnetFiles  = "./networks/HCAS_rect_v%d_pra%d_tau%02d_%dHU.nnet" # File format for .nnet files
modelFiles = "./models/HCAS_rect_v%d_pra%d_tau%02d_%dHU.ckpt" # File format for .nnet files
tbFiles = "./tensorboard/HCAS_rect_v%d_pra%d_tau%02d_%dHU"
##########################

# Custom tensorflow session. Sets up training with either a cpu, gpu, or multiple gpus
def get_session(gpu_ind,gpu_mem_frac=0.45):
    """Create a session that dynamically allocates memory."""
    if gpu_ind[0]>-1:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = ",".join(np.char.mod('%d', gpu_ind))
        config = tf.ConfigProto(device_count = {'GPU': len(gpu_ind)})
        config.gpu_options.per_process_gpu_memory_fraction = gpu_mem_frac
        session = tf.Session(config=config)
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        session = tf.Session()
    return session

# Function to compute network accuracy given. However, this does not 
# add the online costs that were used to generate the correct minimum cost index, so
# these accuracies are only an estimate
def custAcc(y_true,y_pred):
    maxesPred = tf.argmax(y_pred,axis=1)
    inds = tf.argmax(y_true,axis=1)
    diff = tf.cast(tf.abs(inds-maxesPred),dtype='float32')
    ones = tf.ones_like(diff,dtype='float32')
    zeros= tf.zeros_like(diff,dtype='float32')
    l = tf.where(diff<0.5,ones,zeros)
    return tf.reduce_mean(l)


# Write the .nnet files
def writeNNet(params,fileName):      
    #Open the file we wish to write
    with open(fileName,'w') as f2:

        #####################
        # First, we write the header lines:
        # The first line written is just a line of text
        # The second line gives the four values:
        #     Number of hidden layers in the networks
        #     Number of inputs to the networks
        #     Number of outputs from the networks
        #     Maximum size of any hidden layer
        # The third line gives the sizes of each layer, including the input and output layers
        # The fourth line specifies if the network is "symmetric", in which the network was only 
        #     trained on half of the state space and then the other half is mirrored. This
        #     option was explored but not fruitfully, so this value is just set to false, 0.
        # The fifth line specifies the minimum values each input can take.
        # The sixth line specifies the maximum values each input can take.
        #     Inputs passed to the network are truncated to be between this range, since neural
        #     networks are not good at extrapolation
        # The seventh line gives the mean value of each input and of the outputs.
        # The eighth line gives the range of each input and of the outputs.
        #     These two lines are used to map raw inputs to the 0 mean, 1 range of the inputs and outputs
        # used during training
        # The ninth line begins the network weights and biases.
        ####################
        f2.write("// Neural Network File Format by Kyle Julian, Stanford 2016\n")

        #Extract the necessary information and write the header information
        keys = sorted(params.keys())
        keysW = keys[:int(len(keys)/2)]
        keysb = keys[int(len(keys)/2):]
        numLayers = len(keysW)
        inputSize = params[keysW[0]].shape[0]
        outputSize = params[keysb[-1]].shape[0]
        maxLayerSize = inputSize;
        for key in keysW:
            if params[key].shape[1]>maxLayerSize :
                maxLayerSize = params[key].shape[0]

        str = "%d,%d,%d,%d,\n" % (numLayers,inputSize,outputSize,maxLayerSize)
        f2.write(str)
        str = "%d," % inputSize
        f2.write(str)
        for key in keysW:
            str = "%d," % params[key].shape[1]
            f2.write(str)
        f2.write("\n")

        #Write Min, Max, Mean, and Range of each of the inputs on outputs for normalization
        f2.write("0,\n") # COC cost
        f2.write(min_inputs   + "\n") #Minimum Input Values
        f2.write(max_inputs  + "\n") #Maximum Input Values                
        f2.write(means  + "\n") #Means of inputs for normalizations
        f2.write(ranges + "\n") #Ranges of inputs for

        ##################
        # Write weights and biases of neural network
        # First, the weights from the input layer to the first hidden layer are written
        # Then, the biases of the first hidden layer are written
        # The pattern is repeated by next writing the weights from the first hidden layer to the second hidden layer,
        # followed by the biases of the second hidden layer.
        ##################
        for ii in range(len(keysW)):
            data = np.array(params[keysW[ii]]).T
            for i in range(len(data)):
                for j in range(int(np.size(data)/len(data))):
                    str = ""
                    if int(np.size(data)/len(data))==1:
                        str = "%.5e," % data[i] #Five digits written. More can be used, but that requires more more space.
                    else:
                        str = "%.5e," % data[i][j]
                    f2.write(str)
                f2.write("\n")
            data = np.array(params[keysb[ii]]).T
            for i in range(len(data)):
                for j in range(int(np.size(data)/len(data))):
                    str = ""
                    if int(np.size(data)/len(data))==1:
                        str = "%.5e," % data[i] #Five digits written. More can be used, but that requires more more space.
                    else:
                        str = "%.5e," % data[i][j]
                    f2.write(str)
                f2.write("\n")



# Train the Tensorflow model given parameters                
def run_model(session, predict, loss, Xd, yd,writer,vd,saveFile,saveTF,
              epochs=1, batch_size=64,
              training=None, save_every=1, test_size=1000):
    
    # have tensorflow compute accuracy
    accuracy = custAcc(y,predict)
    tf.summary.scalar('accuracy',accuracy)
    tf.summary.scalar('mean_loss',loss)
    merged_summary = tf.summary.merge_all()
    saver = tf.train.Saver()
        
    
    for e in range(epochs):
        # keep track of losses and accuracy
        correct = 0
        
        # shuffle indicies
        train_indices = np.arange(Xd.shape[0])
        np.random.shuffle(train_indices)
        
        # make sure we iterate over the dataset once
        for i in range(int(math.ceil(Xd.shape[0]/batch_size))):
            
            # generate indicies for the batch
            start_idx = (i*batch_size)%Xd.shape[0]
            idx = train_indices[start_idx:np.minimum(start_idx+batch_size,Xd.shape[0])]
            
            # create a feed dictionary for this batch
            feed_dict = {X: Xd[idx,:],
                         y: yd[idx,:]
                        }
            
            # have tensorflow perform a training step
            session.run([training],feed_dict=feed_dict)

        test_inds = np.random.choice(Xd.shape[0],test_size, replace=False)
        # Save to to summary folder, can view data with Tensorboard
        feed_dict = {X: Xd[test_inds,:],
                     y: yd[test_inds,:]
                    }
        s = sess.run(merged_summary, feed_dict=feed_dict)
        writer.add_summary(s,e+1)
            
        # Save graph, write .nnet file
        params = sess.run(vd)
        writeNNet(params,saveFile)
        saver.save(sess, saveTF)
            
        # Save model at specified intervals
        if (e+1) % save_every == 0:
            params = sess.run(vd)
            writeNNet(params,saveFile[:-5] + "_%03d.nnet"%(e+1))
            saver.save(sess, saveTF[:-5] + "_%03d.ckpt"%(e+1))



# The previous RA should be given as a command line input
if len(sys.argv) > 2:
    pra = int(sys.argv[1])
    tau = int(sys.argv[2])
    gpu = -1
    if len(sys.argv)>3:
        gpu = int(sys.argv[3])

    print("Loading Data for HCAS, pra %02d, Network Version %d" % (pra, ver))
    f       = h5py.File(trainingDataFiles % (table_ver,pra,tau),'r')
    X_train = np.array(f['X'])
    Q       = np.array(f['y'])
    means = np.array(f['means'])
    ranges=np.array(f['ranges'])
    mins = np.array(f['min_inputs'])
    maxes = np.array(f['max_inputs'])


    min_inputs      = ",".join(np.char.mod('%f', mins))
    max_inputs     = ",".join(np.char.mod('%f', maxes))
    means      = ",".join(np.char.mod('%f', means))
    ranges     = ",".join(np.char.mod('%f', ranges))
    
    #goodInds = np.where(X_train[:,2]==0.0)[0]
    #X_train = X_train[goodInds,:]
    #X_train = X_train[:,[0,1,3]]
    #Q = Q[goodInds,:]
    #means = means[[0,1,3,4]]
    #ranges = ranges[[0,1,3,4]]
    #min_inputs = min_inputs[[0,1,3]]
    #max_inputs = max_inputs[[0,1,3]]
                  
    
    N,numInputs = X_train.shape
                       
    N,numOut = Q.shape
    print("Setting up Model")

    # Asymmetric loss function
    lossFactor = 40.0
    def asymMSE(y_true, y_pred):
        d = y_true-y_pred
        maxes = tf.argmax(y_true,axis=1)
        maxes_onehot = tf.one_hot(maxes,numOut)
        others_onehot = maxes_onehot-1
        d_opt = d*maxes_onehot 
        d_sub = d*others_onehot
        a = lossFactor*(numOut-1)*(tf.square(d_opt)+tf.abs(d_opt))
        b = tf.square(d_opt)
        c = lossFactor*(tf.square(d_sub)+tf.abs(d_sub))
        d = tf.square(d_sub)
        loss = tf.where(d_sub>0,c,d) + tf.where(d_opt>0,a,b)
        return tf.reduce_mean(loss)

    # Define model architecture
    X = tf.placeholder(tf.float32, [None, numInputs])
    y = tf.placeholder(tf.float32, [None, numOut])
    layer_sizes = [hu, hu, hu, hu, hu]
    layers = np.concatenate(([numInputs],layer_sizes,[numOut]))
    vd = {}
    inputs = X
    for i, (inLayer, outLayer) in enumerate(zip(layers[:-1],layers[1:])):
        vd["W" + str(i)] = tf.get_variable("W"+str(i), shape=[inLayer, outLayer])
        vd["b" + str(i)] = tf.get_variable("b"+str(i), shape=[outLayer])
        if i < len(layers)-2:
            inputs = tf.nn.relu(tf.matmul(inputs,vd["W" + str(i)]) + vd["b" + str(i)])
        else:
            # Don't use ReLU on output layer!
            inputs = tf.matmul(inputs,vd["W" + str(i)]) + vd["b" + str(i)]
    y_out = inputs

    # define our loss
    mean_loss = asymMSE(y,y_out)

    # define our optimizer
    optimizer = tf.train.AdamOptimizer(3e-4) # select optimizer and set learning rate
    train_step = optimizer.minimize(mean_loss)

    # Initialize session
    sess = get_session([gpu], 0.45)
    sess.run(tf.global_variables_initializer())

    tb_file = tbFiles % (ver, pra, tau, hu)

    writer = tf.summary.FileWriter(tb_file)
    writer.add_graph(sess.graph)

    nnet_file = nnetFiles % (ver, pra, tau, hu)
    model_file = modelFiles % (ver, pra, tau, hu)
    
    run_model(sess,y_out,mean_loss,X_train,Q,writer,vd,nnet_file,model_file,totalEpochs,2**9,train_step,saveEvery,1000)

    
    ## Train and write nnet files
    #epoch= saveEvery
    #while epoch <= totalEpochs:
    
    #    model.fit(X_train, Q, nb_epoch=saveEvery, batch_size=2**9,shuffle=True)
    #    saveFile = nnetFiles % (pra, ver,hu,epoch)
    #    saveNNet(model,saveFile,means,ranges,min_inputs,max_inputs)
    #    epoch += saveEvery
