__author__ = 'qiaoqiao'
from math import *
from math import log
from math import exp
from numpy import log1p
import copy
import numpy
from numpy import genfromtxt
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt

# function to calculate log of sums
def sumLogProb(a,b):
    if a>b:
        return a+log1p(exp(b-a))
    else:
        return b+log1p(exp(a-b))


def emission(k,y,mu,theta):
    p=[]

    temp1=1-0.5*(1+erf((0.4-mu[0])/sqrt(2.0*theta[0])))
    p.append(((1.0/(sqrt(2*pi)*theta[0]))*exp(-0.5*((y-mu[0])**2/theta[0])))/temp1)

    temp2=0.5*(1+erf((-0.4-mu[1])/sqrt(2.0*theta[1])))
    p.append(((1.0/(sqrt(2*pi)*theta[1]))*exp(-0.5*((y-mu[1])**2/theta[1])))/temp2)

    p.append(1.0/(sqrt(2*pi)*theta[2])*exp(-0.5*((y-mu[2])**2/theta[2])))
    if p[k]/sum(p)<=10**(-100):
        return 10**(-100)
    else:
        return p[k]/sum(p)


# forward algorithm
# emission are three distribution prob with 6 parameters
def forward(mu,theta,seq,P,a0):
    #a0=[1/4,1/4,1/2] #prob to 1,2,3 --- h,l,n
    f=[]
    # the first row is for high GC content, and the second row is for low GC content
    t1=emission(0,seq[0],mu,theta)*a0[0];
    t2=emission(1,seq[0],mu,theta)*a0[1];
    t3=emission(2,seq[0],mu,theta)*a0[2];
    f.append([log(t1),log(t2),log(t3)])

    for i in range(1,len(seq)):
        y=seq[i]

        t1=sumLogProb(sumLogProb(f[i-1][0]+log(P[0,0]),f[i-1][1]+log(P[1,0])),f[i-1][2]+log(P[2,0]))+log(emission(0,y,mu,theta))
        t2=sumLogProb(sumLogProb(f[i-1][0]+log(P[0,1]),f[i-1][1]+log(P[1,1])),f[i-1][2]+log(P[2,1]))+log(emission(1,y,mu,theta))
        t3=sumLogProb(sumLogProb(f[i-1][0]+log(P[0,2]),f[i-1][1]+log(P[1,2])),f[i-1][2]+log(P[2,2]))+log(emission(2,y,mu,theta))


        f.append([t1,t2,t3])
    return f


# backward algorithm
def backward(mu,theta,seq,P):

    #initial V1
    b=numpy.zeros((3,len(seq)))
    b[:,len(seq)-1]=[0,0,0]
    # the first row is for high GC content, and the second row is for low GC content
    for i in range(len(seq)-2,-1,-1):
        y=seq[i+1]
        t1=sumLogProb(sumLogProb((b[0,i+1]+log(P[0,0])+log(emission(0,y,mu,theta))),(b[1,i+1]+log(P[0,1])+log(emission(1,y,mu,theta)))),(b[2,i+1]+log(P[0,2])+log(emission(2,y,mu,theta))))
        t2=sumLogProb(sumLogProb((b[0,i+1]+log(P[1,0])+log(emission(0,y,mu,theta))),(b[1,i+1]+log(P[1,1])+log(emission(1,y,mu,theta)))),(b[2,i+1]+log(P[1,2])+log(emission(2,y,mu,theta))))
        t3=sumLogProb(sumLogProb((b[0,i+1]+log(P[2,0])+log(emission(0,y,mu,theta))),(b[1,i+1]+log(P[2,1])+log(emission(1,y,mu,theta)))),(b[2,i+1]+log(P[2,2])+log(emission(2,y,mu,theta))))
        b[:,i]=[t1,t2,t3]
    return b


# P is the transition probability matrix
# theta,mu -- vector with theta1, theta2 and theta3
#seq -- matrix each column is a sequence
def EM(trans_prob,theta,mu,seq):
    times=0
    loglikelihood=[]
    (a,b)=seq.shape
    maxlog=float("-inf")
    index=0
    # stop when two successive iterations have difference less than 0.01
    while times<=1 or abs(loglikelihood[times-1]-loglikelihood[times-2])>=0.01:
        # initialization
        P=numpy.zeros((3,3))
        exp_mu=numpy.zeros((2,2))
        exp_theta=numpy.zeros((3,2))
        exp_pi=numpy.zeros((4,1))
        like_seq=[]
        pi=[1.0/4,1.0/4,1.0/2]
        for m in range(b):
            # build emission matrix from theta_h and theta_l
            # calculate forward and backcard probability
            F=forward(mu,theta,seq[:,m],trans_prob,pi)
            B=backward(mu,theta,seq[:,m],trans_prob)
            # calculate the log likelihood
            py=sumLogProb(sumLogProb(F[len(F)-1][0],F[len(F)-1][1]),F[len(F)-1][2])
            like_seq.append(py)
            gamma=[]
            # calculate expectation for transition matrix
            P[0,0]+=sum([exp(F[i][0]+trans_prob[0][0]+log(emission(0,seq[i+1,m],mu,theta))+B[0,i+1]-py) for i in range(0,a-1)])
            P[0,1]+=sum([exp(F[i][0]+trans_prob[0][1]+log(emission(1,seq[i+1,m],mu,theta))+B[1,i+1]-py) for i in range(0,a-1)])
            P[0,2]+=sum([exp(F[i][0]+trans_prob[0][2]+log(emission(2,seq[i+1,m],mu,theta))+B[2,i+1]-py) for i in range(0,a-1)])
            P[1,0]+=sum([exp(F[i][0]+trans_prob[1][0]+log(emission(0,seq[i+1,m],mu,theta))+B[0,i+1]-py) for i in range(0,a-1)])
            P[1,1]+=sum([exp(F[i][0]+trans_prob[1][1]+log(emission(1,seq[i+1,m],mu,theta))+B[1,i+1]-py) for i in range(0,a-1)])
            P[1,2]+=sum([exp(F[i][0]+trans_prob[1][2]+log(emission(2,seq[i+1,m],mu,theta))+B[2,i+1]-py) for i in range(0,a-1)])
            P[2,0]+=sum([exp(F[i][0]+trans_prob[2][0]+log(emission(0,seq[i+1,m],mu,theta))+B[0,i+1]-py) for i in range(0,a-1)])
            P[2,1]+=sum([exp(F[i][0]+trans_prob[2][1]+log(emission(1,seq[i+1,m],mu,theta))+B[1,i+1]-py) for i in range(0,a-1)])
            P[2,2]+=sum([exp(F[i][0]+trans_prob[2][2]+log(emission(2,seq[i+1,m],mu,theta))+B[2,i+1]-py) for i in range(0,a-1)])
            exp_pi[0]+=exp(F[0][0])
            exp_pi[1]+=exp(F[0][1])
            exp_pi[2]+=exp(F[0][2])
            exp_pi[3]+=exp(F[0][0])+exp(F[0][1])+exp(F[0][2])
            # calculate expectation for emission matrix
            for j in range(3):
                gamma.append([exp(F[i][j])*exp(B[j,i])/sum([exp(F[i][k])*exp(B[k,i]) for k in range(3)]) for i in range(a)])
            # update mu, theta_h and theta_l according to expectation
            exp_mu[0,0]+=sum([gamma[0][i]*seq[i,m] for i in range(a)])
            exp_mu[0,1]+=sum(gamma[0])
            exp_mu[1,0]+=sum([gamma[1][i]*seq[i,m] for i in range(a)])
            exp_mu[1,1]+=sum(gamma[1])

            exp_theta[0,0]+=sum([gamma[0][i]*(seq[i,m]-mu[0])**2 for i in range(a)])
            exp_theta[0,1]+=sum(gamma[0])
            exp_theta[1,0]+=sum([gamma[1][i]*(seq[i,m]-mu[1])**2 for i in range(a)])
            exp_theta[1,1]+=sum(gamma[1])
            exp_theta[2,0]+=sum([gamma[2][i]*(seq[i,m]-mu[2])**2 for i in range(a)])
            exp_theta[2,1]+=sum(gamma[2])

            # update emission matrix
        mu[0]=exp_mu[0,0]/exp_mu[0,1]
        mu[1]=exp_mu[1,0]/exp_mu[1,1]
        for k in range(3):
            trans_prob[k,0]=P[k,0]/sum(P[k,:])
            trans_prob[k,1]=P[k,1]/sum(P[k,:])
            trans_prob[k,2]=P[k,2]/sum(P[k,:])
            pi[k]=exp_pi[k]/exp_pi[3]
            #mu[k]=exp_mu[k,0]/exp_mu[k,1]
            theta[k]=exp_theta[k,0]/exp_theta[k,1]

        if sum(like_seq)>=maxlog:
            index=times
            max_mu=copy.deepcopy(mu)
            max_theta=copy.deepcopy(theta)
            max_trans=copy.deepcopy(trans_prob)
            maxlog=sum(like_seq)
            max_pi=pi

        loglikelihood.append(sum(like_seq))

        times+=1
    return max_mu,max_theta,max_trans,max_pi,loglikelihood

my_data = genfromtxt('newchr1.csv', delimiter=',')
(a,b)=my_data.shape
print a,b
seq=my_data[0:a-1,1:21]
print seq.var()

for m in range(20):
    seq[:,m]=seq[:,m]-seq.mean()
plt.hist(seq)
plt.show()

test=[[0.7, 0.03,0.27],[0.08,0.65,0.27],[0.1,0.1,0.8]]
X=numpy.array(test)

[mymu,mytheta,mytrans_prob,mypi,myloglikelihood]=EM(X,[1.0,1.0,1.0],[1,-1,0],seq)
print mymu
print mytheta
print mytrans_prob
print myloglikelihood


def posterior(mu,theta,trans_prob,pi,seq):
    prob=numpy.zeros((3,len(seq)))
    F=forward(mu,theta,seq,trans_prob,pi)
    B=backward(mu,theta,seq,trans_prob)
    py=sumLogProb(sumLogProb(F[len(F)-1][0],F[len(F)-1][1]),F[len(F)-1][2])
    for i in range(len(seq)):
        for j in range(3):
            prob[j,i]=sumLogProb(F[i][j],B[j,i])-py
    return prob

post_prob=posterior(mymu,mytheta,mytrans_prob,mypi,my_data[0:a-1,22])

print post_prob