__author__ = 'qiaoqiao'

from math import log
import math
import numpy as np
import matplotlib.pyplot as plt


# read in sequence matrix
def openfile(filename):
    with open(filename,'r') as data:
        lines = data.read().splitlines()
        i = 0
        array=[]
        while i < len(lines):
            line = lines[i]
            if line.startswith(">"):
                seq = ""
                while i+1 <len(lines) and (not lines[i+1].startswith(">")):
                    seq += lines[i+1]
                    i += 1
                # store sequence in a two-dimensional list
                array.append([str(x) for x in seq])
            i += 1
    return array


# calculate frequency matrix for motif
def cal_freq(X,start,K,ex_seq):
    freq=np.zeros((4,K))

    for j in range(K):
            for i in range(len(X)):
                # skip the current sequence
                if i==ex_seq:
                    continue
                if X[i][int(start[i])+j]=='A':
                    freq[0,j]+=1
                if X[i][int(start[i])+j]=='C':
                    freq[1,j]+=1
                if X[i][int(start[i])+j]=='G':
                    freq[2,j]+=1
                if X[i][int(start[i])+j]=='T':
                    freq[3,j]+=1
    return freq


# gibbs sampling algorithm
def gibbs(K,seq,epsilon,times):
    index={'A':0,'C':1,'G':2,'T':3}
    t=len(seq)
    # initial start positions for each motif

    loglike=[]
    indicator=0
    max_log=-np.inf
    max_start = []
    # running some times for a fix number of iterations
    while indicator <10:
        print 'Indicator: '+str(indicator)
        s=np.zeros((K,1))
        for i in range(t):
            test=np.random.choice((len(seq[i])-K+1), 1)
            s[i]=test[0]
        new_start=s
        # iterations for sampling
        for n in range(times):
            # calculate weight matrix for each sequence
            for i in range(t):
                # calculate frequency matrix for new set of start positions without sequence i
                freq_matrix=cal_freq(seq,new_start,K,i)

                # calculate weight matrix from frequency matrix and exclude current sequence
                weight_matrix=(freq_matrix+epsilon)*1.0/(t-1+4*epsilon)
                weight_start=[]
                # for a given start position r, calculate the weight
                for r in range(len(seq[i])-K+1):
                    temp=np.prod([weight_matrix[index[seq[i][r+l]],l] for l in range(K)])
                    weight_start.append(temp)
                # normalize the weight to probability
                prob=[weight_start[j]/sum(weight_start) for j in range(len(weight_start))]
                # sample new start position in proportion to weights
                test=np.random.choice(len(prob), 1, p=prob)
                new_start[i]=test[0]
                
            # calculate the log likelihood for new site
            pi=cal_freq(seq,new_start,K,-1)
            pi=(pi+epsilon)/(t+4*epsilon)
            part=0
            for i in range(t):
                for j in range(K):
                    part+=log(pi[index[seq[i][int(new_start[i])+j]],j])
                part +=(len(seq[i])-K)*log(0.25)

            # update maximum log likelihood and corresponding weights and start positions
            if part>max_log:
                print 'Updated'
                max_log=part
                max_start=new_start
        
        loglike.append(max_log)
        indicator +=1
        #print indicator
    return max_start,max_log

k=10
array1=openfile('motif1.fa')
array2=openfile('motif2.fa')
# run the algorithm for motif1.fa to find motif
start1,likelihood1=gibbs(k,array1,0.25,500)

#motif=np.argmax(result1,axis=0)
#map={0:'A',1:'C',2:'G',3:'T'}
#print 'The motif with maximum likelihood:'
#print ''.join([map[motif[i]] for i in range(len(motif))])
print 'The motif in sequences:'
for i in range(len(start1)):
    print ''.join(array1[i][int(start1[i]):(int(start1[i])+k)])

print start1
print likelihood1
#plt.figure(1)
#plt.plot(likelihood1)

# start2,likelihood2=gibbs(k,array2,0.25,500)

# #motif=np.argmax(result1,axis=0)
# #map={0:'A',1:'C',2:'G',3:'T'}
# #print 'The motif with maximum likelihood:'
# #print ''.join([map[motif[i]] for i in range(len(motif))])
# print 'The motif in sequences:'
# for i in range(len(start2)):
#     print ''.join(array2[i][int(start2[i]):(int(start2[i])+k)])

# print start2
# print likelihood2


# run the algorithm for motif2.fa to find motif
# result2,start2,likelihood2=gibbs(k,array2,0.25,500)
# print result2
# motif=np.argmax(result2,axis=0)

# print 'The motif with maximum likelihood:'
# print ''.join([map[motif[i]] for i in range(len(motif))])
# print 'The motif in sequences:'
# for i in range(len(start2)):
#     print ''.join(array2[i][int(start2[i]):(int(start2[i])+k)])

# print start2
# print likelihood2
# plt.figure(2)
# plt.plot(likelihood2)
# plt.show()
