#!/usr/bin/env python

import os, sys, math, copy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#######################					
# Event class
#######################
class Event:

	def __init__(self):
		self.decayTheta = -99.
		self.decayPhi = -99.
		self.decayR = -99.
		self.totHits = 0
		self.numHitPMTs = 0
		self.hitPMTs = []
		
	def SetSummary(self,decayTheta,decayPhi,decayR,totHits,numHitPMTs):
		self.decayTheta = float(copy.deepcopy(decayTheta))
		self.decayPhi = float(copy.deepcopy(decayPhi))
		self.decayR = float(copy.deepcopy(decayR))
		self.totHits = int(copy.deepcopy(totHits))
		self.numHitPMTs = int(copy.deepcopy(numHitPMTs))
		
	def AddHitPMT(self,hit):
		self.hitPMTs.append(copy.deepcopy(hit))


########################
# Main program
########################

events = []
# read file with data
with open('data.txt', 'r') as f:
	for line in f:
		event = Event()
		splitline = line.split('\t')
		event.SetSummary(splitline[0],splitline[1],splitline[2],splitline[3],splitline[4])
		for i in range(5,len(splitline)-1):
			event.AddHitPMT(map(int,splitline[i].split(' ')))
		events.append(event)
		
f.closed

# prepare plots by generating arrays of quantities to plot
numPMTs = []
for event in events:
	numPMTs.append(float(event.numHitPMTs))
		
# do the plots
plt.figure(figsize=(12,6))
gs = gridspec.GridSpec(1, 2)

# subfigure 1 (histogram of number of hit PMTs)
plt.subplot(gs[0])
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.xlabel('number of hit PMTs',fontsize=16)
plt.ylabel('#', fontsize=16)

h1 = plt.hist(numPMTs,10,facecolor='b',label='number of hit PMTs')

# finish plotting
plt.show(block=False)
raw_input("Press any Key to continue")
plt.close()