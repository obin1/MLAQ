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

A = np.matrix('0, 1,-1, 0, 0, 0, 0, 0, 0, 0;\
           1, 0,-1, 0, 0, 0,-1, 0, 0, 0;\
          -1, 0, 1, 0, 0, 0, 1,-1, 0, 0;\
           0, 0, 0,-1,-1,-1, 0, 0, 0, 0;\
           0, 0, 0, 2, 0, 1,-1, 0, 0, 1;\
           0, 0, 0, 0, 0, 0, 0, 0,-1,-1;\
           0, 0, 0, 0, 0,-1, 1,-1, 2,-1;\
           1,-1, 0, 0, 0, 0, 0, 0, 0, 0;\
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0;\
           0, 0, 0, 1, 1, 1, 0, 0, 0, 0;\
           0, 0, 0, 0, 1, 0, 0, 0, 0, 0')
Ared = np.matrix('0, 1,-1, 0, 0, 0, 0, 0, 0, 0;\
           1, 0,-1, 0, 0, 0,-1, 0, 0, 0;\
          -1, 0, 1, 0, 0, 0, 1,-1, 0, 0;\
           0, 0, 0,-1,-1,-1, 0, 0, 0, 0;\
           0, 0, 0, 2, 0, 1,-1, 0, 0, 1;\
           0, 0, 0, 0, 0, 0, 0, 0,-1,-1')

#D = pd.read_csv('/Users/obin/Dropbox/ML Share/KleemanFortran/base/delC.txt',header=None,delim_whitespace=True).values
#S = pd.read_csv('/Users/obin/Dropbox/ML Share/KleemanFortran/base/S.txt',header=None,delim_whitespace=True).values
#C = pd.read_csv('/Users/obin/Dropbox/ML Share/KleemanFortran/base/C.txt',header=None,delim_whitespace=True).values
#J = pd.read_csv('/Users/obin/Dropbox/ML Share/KleemanFortran/base/J.txt',header=None,delim_whitespace=True).values
#day = pd.read_csv('/Users/obin/Dropbox/ML Share/KleemanFortran/base/day.txt',header=None,delim_whitespace=True).values
#time = pd.read_csv('/Users/obin/Dropbox/ML Share/KleemanFortran/base/time.txt',header=None,delim_whitespace=True).values

#D = pd.read_csv('/Users/obin/Dropbox/ML Share/MLAQ/delC.txt',header=None,delim_whitespace=True).values
S = pd.read_csv('/Users/obin/Dropbox/ML Share/MLAQ/S.txt',header=None,delim_whitespace=True).values
C = pd.read_csv('/Users/obin/Dropbox/ML Share/MLAQ/C.txt',header=None,delim_whitespace=True).values
J = pd.read_csv('/Users/obin/Dropbox/ML Share/MLAQ/J.txt',header=None,delim_whitespace=True).values
day = pd.read_csv('/Users/obin/Dropbox/ML Share/MLAQ/day.txt',header=None,delim_whitespace=True).values
time = pd.read_csv('/Users/obin/Dropbox/ML Share/MLAQ/time.txt',header=None,delim_whitespace=True).values


D = np.delete(np.diff(C,axis = 0), list(range(240, C.shape[0], 241)), axis=0)

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
    if np.allclose(D[i,:],A@S[i,:], rtol = 1e-8):
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
trainD = D[np.logical_or(day[:,0] <= 500,day[:,0] > 600),:]
trainS = S[np.logical_or(day[:,0] <= 500,day[:,0] > 600),:]
trainC = startC[np.logical_or(day[:,0] <= 500,day[:,0] > 600),:]
trainJ = startJ[np.logical_or(day[:,0] <= 500,day[:,0] > 600)]
traintime = time[np.logical_or(day[:,0] <= 500,day[:,0] > 600)]
trainday = day[np.logical_or(day[:,0] <= 500,day[:,0] > 600)]
x_train = np.hstack([trainC[:,0:6],trainJ])
scaler = MinMaxScaler()
x_train = scaler.fit_transform(x_train)

testD = D[np.logical_and(day[:,0] > 500,day[:,0] <= 600),:]
testS = S[np.logical_and(day[:,0] > 500,day[:,0] <= 600),:]
testC = startC[np.logical_and(day[:,0] > 500,day[:,0] <= 600),:]
testJ = startJ[np.logical_and(day[:,0] > 500,day[:,0] <= 600)]
testtime = time[np.logical_and(day[:,0] > 500,day[:,0] <= 600)]
testday = day[np.logical_and(day[:,0] > 500,day[:,0] <= 600)]
x_test = np.hstack([testC[:,0:6], testJ])
x_test = scaler.transform(x_test)

# Train NN on S values
o3noxrxns = [0, 1, 2, 6, 7]
singlerxns = [3, 4, 5, 8, 9]
model = [None]*6
model[0] = Sequential()
model[0].add(Dense(20, activation='relu'))
model[0].add(Dropout(0.2))
model[0].add(Dense(40, activation='relu'))
model[0].add(Dropout(0.2))
model[0].add(Dense(40, activation='relu'))
model[0].add(Dropout(0.2))
model[0].add(Dense(20, activation='relu'))
model[0].add(Dropout(0.2))
model[0].add(Dense(5, activation='sigmoid'))
model[0].compile(loss='mean_squared_error',
                  optimizer='adam',
                  metrics=['mean_absolute_error', 'mean_squared_error'])
#model[0].summary()
#    model[r].fit(x_train, trainS[:,2], epochs=10, batch_size=32)
model[0].fit(x_train, trainS[:,o3noxrxns]*1000, epochs=10, batch_size=8)

for r in range(0,np.shape(singlerxns)[0]):
    print('training NN for reaction ', singlerxns[r]+1)
    model[r+1] = Sequential()
    model[r+1].add(Dense(10, activation='relu', input_dim=7))
    model[r+1].add(Dropout(0.2))
    model[r+1].add(Dense(20, activation='relu'))
    model[r+1].add(Dropout(0.2))
    model[r+1].add(Dense(1, activation='sigmoid'))
    model[r+1].compile(loss='mean_squared_error',
                  optimizer='adam',
                  metrics=['mean_absolute_error', 'mean_squared_error'])
#    model[r+1].summary()
#    model[r].fit(x_train, trainS[:,2], epochs=10, batch_size=32)
    model[r+1].fit(x_train, trainS[:,singlerxns[r]]*1000, epochs=5, batch_size=8)

# Predict S 
Spredict = np.zeros(np.shape(testS))
Dpredict = np.zeros(np.shape(testD))
Spredict[:,o3noxrxns] = model[0].predict(x_test)/1000
for r in range(0,np.shape(singlerxns)[0]):
    Spredict[:,singlerxns[r]:singlerxns[r]+1] = model[r+1].predict(x_test)/1000
  
# Make delC predictions from concentration over 6 minutes

for i in range (0,np.shape(Spredict)[0]):
    Dpredict[i,:] = A@Spredict[i,:]
Cpredict = testC + Dpredict

# Given a starting concentration, predict over a day
    # Note: this should be made its own function
chosenday = 505 #505 (507, 567 nice too)
locs = np.where(testday == chosenday)[0]
#conc = np.expand_dims(rawC[241*(chosenday-1),0:6], axis = 0)
conc = np.expand_dims(testC[locs[0],:], axis = 0)
daypredict = np.zeros([240,11])
sguess = np.zeros([10,1])
for t in range(0,240):
#    if (t%1 == 0):
##        print(t)
#        conc = np.expand_dims(rawC[241*(chosenday-1)+t,:], axis = 0)
    cond = np.hstack([conc[:,0:6], np.expand_dims(rawJ[241*(chosenday-1)+t], axis =0)])
    cond = scaler.transform(cond)
    sguess[o3noxrxns] = np.transpose(model[0].predict(cond)/1000)
    for r in range(0,np.shape(singlerxns)[0]):
        sguess[singlerxns[r]] = model[r+1].predict(cond)/1000
#    sguess[o3noxrxns,0] = rawS[240*(chosenday-1)+t,o3noxrxns]
    ## Check to see if C(t+1) = C(t) + A*S
#    sguess = rawS[240*(chosenday-1)+t,0:10] 
#    if sguess[0] > 7776:
#        print(t)
#        sguess[0] = 0
    delC = A@sguess
#    delC = rawD[240*(chosenday-1)+t,:]
    conc = conc+delC.transpose()
#    conc[conc <0] = 0;
    daypredict[t,:] = conc


# Plotting time
compound = ['O3','NO','NO2','HCHO','HO2','HO2H']
for c in range(0,6):
    plt.figure(c)
    plt.plot(testtime[locs[1:]]/60,testC[locs[1:],c], 'b-', label = 'Processed Fortran', linewidth =0.8)
    plt.plot((testtime[locs[0:-1]]+6)/60,Cpredict[locs[0:-1],c], 'r-', label = 'NN 6 minute intervals', linewidth =0.8)
    plt.plot(np.arange(testtime[locs[1]]/60,24.09,0.1),daypredict[:,c], 'g-', label = 'NN Full Day', linewidth =0.8)
    plt.plot(np.arange(0.1,24.1,0.1),rawC[241*(chosenday-1):241*(chosenday)-1,c], label = 'Raw fort', linewidth =0.8)
    plt.ylabel('Concentration [ppm]')
    plt.xlabel('Time [hours]')
    plt.legend()
    plt.title(compound[c])
#    plt.savefig(compound[c] + '.png', dpi=1200)

totalCarbon_NN = daypredict[:,3] + daypredict[:,9] # [HCHO] + [CO]
totalNitrogen_NN = daypredict[:,1] + daypredict[:,2] + daypredict[:,8] # [NO] + [NO2] + [HNO3] 
totalCarbon_FT = rawC[241*(chosenday-1):241*(chosenday)-1,3] + rawC[241*(chosenday-1):241*(chosenday)-1,9] # [HCHO] + [CO]
totalNitrogen_FT = rawC[241*(chosenday-1):241*(chosenday)-1,1] + rawC[241*(chosenday-1):241*(chosenday)-1,2] + rawC[241*(chosenday-1):241*(chosenday)-1,8] # [NO] + [NO2] + [HNO3] 

# np.std(totalNitrogen_FT) ~ 4.17e-5 with 1e-3 | ~ 0.000565 with 2e-4 timesteps
# np.std(totalCarbon_FT) ~ 1.75e-5 with 1e-3   | ~ 0.000332 with 2e-4 timesteps
 
plt.figure(6)
plt.plot(np.arange(testtime[locs[1]]/60,24.09,0.1),totalCarbon_NN, label = 'Neural Network')
plt.plot(np.arange(testtime[locs[1]]/60,24.09,0.1), totalCarbon_FT, label = 'Raw Fortran')
plt.ylabel('Total Carbon Concentration [ppm]')
plt.xlabel('Time [hours]')
plt.legend()
plt.title('day ' + str(chosenday))
#plt.savefig('TotalCarbon_Julia.png', dpi=1200)

plt.figure(7)
plt.plot(np.arange(testtime[locs[1]]/60,24.09,0.1),totalNitrogen_NN, label = 'Neural Network')
plt.plot(np.arange(testtime[locs[1]]/60,24.09,0.1), totalNitrogen_FT, label = 'Raw Fortran')
plt.ylabel('Total Nitrogen Concentration [ppm]')
plt.xlabel('Time [hours]')
plt.legend()
plt.title('day ' + str(chosenday))
#plt.savefig('TotalNitrogen_Julia.png', dpi=1200)


    
    
#plt.plot(testtime.iloc[np.where(testday.iloc[:,0] == 301)[0]],testC.iloc[np.where(testday.iloc[:,0] == 301)[0],0])
#plt.plot(testtime[locs[1:-1]]/60,testJ[locs[1:-1]], label = 'act flux')
  