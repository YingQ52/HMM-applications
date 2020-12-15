__author__ = 'qiaoqiao'
import csv
from math import *
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
            line = [x for x in lines[i].split()]
            array.append(line)
            i += 1
        k=0
        while k < len(array[0]):
            temp=[array[j][k] for j in range(len(array))]
            if temp.count("NA")>len(line)/1000.0:
                for row in array:
                    del row[k]
            k+=1
        print len(array)
        print len(array[0])

        result=[]
        for i in range(len(array)):
            if array[i].count("NA") ==0:
                line = [float(x) for x in array[i]]
                result.append(line)
        #         # store sequence in a two-dimensional list

    return result


# chr1=openfile('chr1.d2out')
# chr2=openfile('chr2.d2out')
# chr3=openfile('chr3.d2out')
# chr4=openfile('chr4.d2out')
# chr5=openfile('chr5.d2out')
# chr6=openfile('chr6.d2out')
# chr7=openfile('chr7.d2out')
chr8=openfile('chr8.d2out')
# chr9=openfile('chr9.d2out')
# chr10=openfile('chr10.d2out')
# chr11=openfile('chr11.d2out')
# chr12=openfile('chr12.d2out')
# chr13=openfile('chr13.d2out')
# chr14=openfile('chr14.d2out')
# chr15=openfile('chr15.d2out')
# chr16=openfile('chr16.d2out')
# chr17=openfile('chr17.d2out')
# chr18=openfile('chr18.d2out')
# chr19=openfile('chr19.d2out')
# chr20=openfile('chr20.d2out')
# chr21=openfile('chr21.d2out')
# chr22=openfile('chr22.d2out')
# #
#
# with open("chr1.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr1)
# with open("chr2.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr2)
# with open("chr3.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr3)
# with open("chr4.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr4)
# with open("chr5.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr5)
# with open("chr6.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr6)
# with open("chr7.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr7)
# with open("chr8.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr8)
# with open("chr9.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr9)
# with open("chr10.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr10)
# with open("chr11.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr11)
# with open("chr12.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr12)
# with open("chr13.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr13)
# with open("chr14.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr14)
# with open("chr15.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr15)
# with open("chr16.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr16)
# with open("chr17.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr17)
# with open("chr18.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr18)
# with open("chr19.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr19)
# with open("chr20.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr20)
# with open("chr21.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr21)
# with open("chr22.csv", "wb") as f:
#     writer = csv.writer(f)
#     writer.writerows(chr22)
#
# x = chr1[3]
# >>> y = np.sin(x)
# >>> xvals = np.linspace(0, 2*np.pi, 50)
# >>> yinterp = np.interp(xvals, x, y)