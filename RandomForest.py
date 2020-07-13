#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:17:02 2019

@author: obin
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 23:01:16 2019

@author: obin
"""

#  Species
#         1) O3
#         2) NO
#         3) NO2
#         4) HCHO
#         5) HO2
#         6) HO2H
#         steady state species: .HO, O
#         buildup species: HNO3, CO, H2

#  Reactions                         
#         1) NO2 + HV = NO + O                              
#         2) O + O2 = O3
#         3) O3 + NO = NO2 + O2
#         4) HCHO + HV = 2 HO2. + CO
#         5) HCHO + HV = H2 + CO
#         6) HCHO + HO. = HO2. + CO + H2O
#         7) HO2. + NO = HO. + NO2
#         8) HO. + NO2 = HNO3
#         9) HO2H + HV = 2 HO.
#        10) HO2H + HO. = H2O + HO2.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, LSTM
from tensorflow.keras.optimizers import SGD
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestRegressor

A = np.matrix('0, 1,-1, 0, 0, 0, 0, 0, 0, 0;\
           1, 0,-1, 0, 0, 0,-1, 0, 0, 0;\
          -1, 0, 1, 0, 0, 0, 1,-1, 0, 0;\
           0, 0, 0,-1,-1,-1, 0, 0, 0, 0;\
           0, 0, 0, 2, 0, 1,-1, 0, 0, 1;\
           0, 0, 0, 0, 0, 0, 0, 0,-1,-1')

D = pd.read_csv('delC.txt',header=None,delim_whitespace=True).values
S = pd.read_csv('S.txt',header=None,delim_whitespace=True).values
C = pd.read_csv('C.txt',header=None,delim_whitespace=True).values
J = pd.read_csv('J.txt',header=None,delim_whitespace=True).values
day = pd.read_csv('day.txt',header=None,delim_whitespace=True).values
time = pd.read_csv('time.txt',header=None,delim_whitespace=True).values

rawday = day;
rawJ = J;
rawC = C;
rawS = S;
rawD = D;
rawtime = time;

startC = np.delete(C, list(range(240, C.shape[0], 241)), axis=0)
startJ = np.delete(J, list(range(240, J.shape[0], 241)), axis=0)
day = np.delete(day, list(range(240, day.shape[0], 241)), axis=0)
time = np.delete(time, list(range(240, time.shape[0], 241)), axis=0)

#plt.plot(C.iloc[1:240,4])

D = D[S[:,0]<7770]
startC = startC[S[:,0]<7770]
startJ = startJ[S[:,0]<7770]
day = day[S[:,0]<7770]
time = time[S[:,0]<7770]
S = S[S[:,0]<7770]

D = D[startC[:,0]<0.2,:]
startJ = startJ[startC[:,0]<0.2,:]
day = day[startC[:,0]<0.2,:]
time = time[startC[:,0]<0.2,:]
S = S[startC[:,0]<0.2,:]
startC = startC[startC[:,0]<0.2,:]

# Find good runs (where A*S = delC)
goodrun = [0]
badrun = [0]
for i in range(0,np.shape(D)[0]):
    if np.allclose(D[i,0:6],A@S[i,:], atol = 1e-9):
        goodrun = np.vstack([goodrun, i])
    else: 
        badrun = np.vstack([badrun, i])
            
goodrun = goodrun[1:-1]
badrun = badrun[1:-1]

D = D[goodrun]
startC = startC[goodrun]
startJ = startJ[goodrun]
day = day[goodrun]
time = time[goodrun]
S = S[goodrun]



#numdays = np.shape(D)
#np.shape(trainD)
#
trainD = D[day[:,0] <= 500,0:6]
trainS = S[day[:,0] <= 500,:]
trainC = startC[day[:,0] <= 500,0:6]
trainJ = startJ[day[:,0] <= 500]
traintime = time[day[:,0] <= 500]
trainday = day[day[:,0] <= 500]
x_train = np.hstack([trainC,trainJ])
scaler = MinMaxScaler()
x_train = scaler.fit_transform(x_train)

testD = D[day[:,0] > 500,0:6]
testS = S[day[:,0] > 500,:]
testC = startC[day[:,0] > 500,0:6]
testJ = startJ[day[:,0] > 500,:]
testtime = time[day[:,0] > 500]
testday = day[day[:,0] > 500]
x_test = np.hstack([testC, testJ])
x_test = scaler.transform(x_test)

# Train NN on S values
o3noxrxns = [0, 1, 2, 6, 7]
forest = [None]*6
forest[0] = RandomForestRegressor(n_estimators = 100)
forest[0].fit(x_train, trainS[:,o3noxrxns]*1000)

singlerxns = [3, 4, 5, 8, 9]
for r in range(0,np.shape(singlerxns)[0]):
    forest[r+1] = RandomForestRegressor(n_estimators = 100)
    print('training RForest for reaction ', singlerxns[r]+1)
    forest[r+1].fit(x_train, trainS[:,singlerxns[r]]*1000)

# Predict S 
Spredict = np.zeros(np.shape(testS))
Dpredict = np.zeros(np.shape(testD))
Spredict[:,o3noxrxns] = forest[0].predict(x_test)/1000
for r in range(0,np.shape(singlerxns)[0]):
    Spredict[:,singlerxns[r]] = forest[r+1].predict(x_test)/1000
  
# Make delC predictions from concentration over 6 minutes

for i in range (0,np.shape(Spredict)[0]):
    Dpredict[i,:] = A@Spredict[i,:]
Cpredict = testC + Dpredict

# Given a starting concentration, predict over a day
    # Note: this should be made its own function
chosenday = 505#505
locs = np.where(testday == chosenday)[0]
#conc = np.expand_dims(rawC[241*(chosenday-1),0:6], axis = 0)
conc = np.expand_dims(testC[locs[0],:], axis = 0)
daypredict = np.zeros([240,6])
sguess = np.zeros([10,1])
for t in range(0,240):
#    if (t%1 == 0):
#        conc = np.expand_dims(rawC[241*(chosenday-1)+t,0:6], axis = 0)
    cond = np.hstack([conc, np.expand_dims(rawJ[241*(chosenday-1)+t], axis =0)])
    cond = scaler.transform(cond)
    sguess[o3noxrxns] = np.transpose(forest[0].predict(cond)/1000)
    for r in range(0,np.shape(singlerxns)[0]):
        sguess[singlerxns[r]] = forest[r+1].predict(cond)/1000
#    sguess[o3noxrxns,0] = rawS[240*(chosenday-1)+t,o3noxrxns]
    ## Check to see if C(t+1) = C(t) + A*S
#    sguess = rawS[240*(chosenday-1)+t,0:10] 
#    if sguess[0] > 7776:
#        print(t)
#        sguess[0] = 0
    delC = A@sguess
#    delC = rawD[240*(chosenday-1)+t,0:6]
    conc = conc+delC.transpose()
#    conc[conc <0] = 0;
    daypredict[t,:] = conc


# Plotting time
compound = ['O3','NO','NO2','HCHO','HO2','HO2H']
for c in range(0,6):
    plt.figure(c)
    plt.plot(testtime[locs[1:]]/60,testC[locs[1:],c], 'b-', label = 'Processed Fortran', linewidth =0.8)
    plt.plot((testtime[locs[0:-1]]+6)/60,Cpredict[locs[0:-1],c], 'r-', label = 'NN 6 minute intervals', linewidth =0.8)
    plt.plot(np.arange(testtime[locs[1]]/60,24.1,0.1),daypredict[:,c], 'g-', label = 'NN Full Day', linewidth =0.8)
#    plt.plot(np.arange(0.1,24.1,0.1),rawC[241*(chosenday-1):241*(chosenday)-1,c], label = 'Raw fort')
    plt.ylabel('Concentration [ppm]')
    plt.xlabel('Time [hours]')
    plt.legend()
    plt.title(compound[c])
#    plt.savefig(compound[c] + '.png', dpi=1200)

totalCarbon_NN = daypredict[:,3] + daypredict[:,9] # [HCHO] + [CO]
totalNitrogen_NN = daypredict[:,1] + daypredict[:,2] + daypredict[:,8] # [NO] + [NO2] + [HNO3] 
totalCarbon_FT = rawC[241*(chosenday-1):241*(chosenday)-1,3] + rawC[241*(chosenday-1):241*(chosenday)-1,9] # [HCHO] + [CO]
totalNitrogen_FT = rawC[241*(chosenday-1):241*(chosenday)-1,1] + rawC[241*(chosenday-1):241*(chosenday)-1,2] + rawC[241*(chosenday-1):241*(chosenday)-1,8] # [NO] + [NO2] + [HNO3] 

plt.figure(11)
plt.plot(np.arange(testtime[locs[1]]/60,24.1,0.1),totalCarbon_NN)
plt.plot(np.arange(testtime[locs[1]]/60,24.1,0.1), totalCarbon_FT)
    
#plt.plot(testtime.iloc[np.where(testday.iloc[:,0] == 301)[0]],testC.iloc[np.where(testday.iloc[:,0] == 301)[0],0])
#plt.plot(testtime[locs[1:-1]]/60,testJ[locs[1:-1]], label = 'act flux')
  