__author__ = 'qiaoqiao'
from math import log
from math import exp
from numpy import log1p
import numpy
import matplotlib.pyplot as plt


# function to calculate log of sums
def sumLogProb(a,b):
    if a>b:
        return a+log1p(exp(b-a))
    else:
        return b+log1p(exp(a-b))


# forward_backward algorithm to calculate posterior probability
def for_backward(emission,seqfile,mu):
    # read in sequence data
    with open(seqfile, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            if line == '\n':
                break
            else:
                seq=[s for s in str(line.split()[0])]
    # get the forward probability matrix from forward algorithm
    F=forward(emission,seq,mu)
    # get the backward probability matrix from backward algorithm
    B=backward(emission,seq,mu)
    # get the posterior probability from the termination step of forward algorithm
    py1=sumLogProb(F[len(F)-1][0],F[len(F)-1][1])
    # get the posterior probability from the termination step of backward algorithm
    py2=sumLogProb((B[0,0]+log(0.5)+log(emission[tuple(['h',seq[0]])])),(B[1,0]+log(0.5)+log(emission[tuple(['l',seq[0]])])))
    post=[]
    py=[]
    for i in range(len(seq)):
        temp=sumLogProb(F[i][0]+B[0,i],F[i][1]+B[1,i])
        py.append(temp)
        post.append([F[i][0]+B[0,i]-temp,F[i][1]+B[1,i]-temp])

    return post,py,py1,py2


# forward algorithm
def forward(emission,seq,mu):
    f=[]

    t1=emission[tuple(['h',seq[0]])]*0.5;
    t2=emission[tuple(['l',seq[0]])]*0.5;
    f.append([log(t1),log(t2)])

    for i in range(1,len(seq)):
        y=seq[i]
        t1=sumLogProb(f[i-1][0]+log(1-mu),f[i-1][1]+log(mu))+log(emission[tuple(['h',y])])
        t2=sumLogProb(f[i-1][0]+log(mu),f[i-1][1]+log(1-mu))+log(emission[tuple(['l',y])])
        f.append([t1,t2])
    return f


# backward algorithm
def backward(emission,seq,mu):

    #initial V1
    b=numpy.zeros((2,len(seq)))
    b[:,len(seq)-1]=[0,0]

    for i in range(len(seq)-2,-1,-1):
        y=seq[i+1]
        t1=sumLogProb((b[0,i+1]+log(1-mu)+log(emission[tuple(['h',y])])),(b[1,i+1]+log(mu)+log(emission[tuple(['l',y])])))
        t2=sumLogProb((b[0,i+1]+log(mu)+log(emission[tuple(['h',y])])),(b[1,i+1]+log(1-mu)+log(emission[tuple(['l',y])])))

        b[:,i]=[t1,t2]

    return b


emission={('h','A'):0.13,('h','C'):0.37,('h','G'):0.37,('h','T'):0.13,
          ('l','A'):0.32,('l','C'):0.18,('l','G'):0.18,('l','T'):0.32}

p1,p2,p3,p4=for_backward(emission,'pset3-sequence.fa',0.01)
ph=[exp(p1[i][0]) for i in range(len(p1))]

print 'marginal posterior probability of state h at every position in the sequence:'
print ph
print 'posterior probabilities P(y) for each position:'
print p2
print 'posterior probabilities P(y) from forward algorithm:'
print p3
print 'posterior probabilities P(y) from backward algorithm:'
print p4


# plot the result
GC=[[66, 415], [527, 720], [951, 1000]]
fig, ax = plt.subplots()
ax.plot(list(range(1,1001,1)),ph)
ax.axvspan(GC[0][0],GC[0][1], alpha=0.1, color='red')
ax.axvspan(GC[1][0],GC[1][1], alpha=0.1, color='red')
ax.axvspan(GC[2][0],GC[2][1], alpha=0.1, color='red')

plt.axis([1, 1000, -0.5, 1.5])
plt.ylabel('probability')
plt.xlabel('position in sequence')
plt.title('probability VS sequence position')
plt.savefig('posterior.png')
plt.show()







