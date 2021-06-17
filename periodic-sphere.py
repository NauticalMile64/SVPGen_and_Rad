'''
The MIT License (MIT)
Copyright (c) 2016 Nolan Dyck

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

This code has been run successfully on Ubuntu 14.10 with yadedaily version 1.20.0-72-bbff50f~trusty

More recently tested with yade 2020.01a on Ubuntu 20.04.2 LTS with python 3.8.5

I am currently getting segfault errors when the GlobalStiffnessTimeStepper is active, so I have deactivated it currently. I was using it to dynamically change the time step to speed convergence, but it appears not to be working for this version of yade.
'''

from yade import pack,export,utils
import math,threading
import _pickle as cPickle
import pickle as pkl
import numpy as np

maxBadRun=5	#Maximum number of times the simulation is allowed to restart before exiting with failure
alMaxF=5e-2	#Divergence criteria: simulation stops and restarts if the average penetration force exceeds this value
alAvgF=5e-3	#Divergence criteria: simulation stops and restarts if the average penetration force exceeds this value
dtStart=8e-6	#Starting time step size (may be decreased if the simulation is in a bad state)

porCheckInterval = 5000		#How often the porosity is checked for convergence
checkPointInterval = 10000	#How often the simulation is checkpointed for potential restarts

r=3e-4			#Mean Radius (m)
s=0.5		#Sphere radii are uniformly distributed between r*(1-s) and r*(1+s)
N=100			#Number of spheres
t=0.065		#Surface tension coefficient (N/m)
p=0.8			#Desired porosity
c=0
rMin=1.0e-6	#Spheres with radii below this limit are excluded

utils.readParamsFromTable(
	Row=0,
	Por=p,
	rMean=r,
	rStd=s,
	sig=t,
	numSpheres=N,
	Save=1
)

from yade.params.table import *
tPorosity=Por
REVL=math.pow(5*numSpheres*(float(4)/3*math.pi*rMean**3),float(1)/3) 
O.dt=dtStart

done=False

#Create the pack of spheres
O.materials.append(BubbleMat(density=1e5))
sp=pack.SpherePack()
tnum=sp.makeCloud((0,0,0),(REVL,REVL,REVL),rMean=rMean,rRelFuzz=s,periodic=True,num=numSpheres)
sphVol=sp.relDensity()*REVL**3
tStrain=-10
REVl=0

for sphere in sp:
	if(sphere[1] < rMin):
		sp.remove(sphere)
		numSpheres -= 1

sp.toSimulation() 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb()],verletDist=0.5*r,allowBiggerThanPeriod=True),
	InteractionLoop([Ig2_Sphere_Sphere_ScGeom()],[Ip2_BubbleMat_BubbleMat_BubblePhys()],[Law2_ScGeom_BubblePhys_Bubble(pctMaxForce=0.15,surfaceTension=sig)]),
	NewtonIntegrator(damping=0.01),
	PeriTriaxController(goal=(tStrain,tStrain,tStrain),stressMask=0b000,maxUnbalanced=alMaxF, globUpdate=1,doneHook='reachedTargetStrain()',label='triax',dynCell=True,mass=1e8),
	GlobalStiffnessTimeStepper(active=False),#,defaultDt=dtStart,targetDt=dtStart,timestepSafetyCoefficient=1,label='stepper', timeStepUpdateInterval=checkPointInterval*100,maxDt=dtStart,densityScaling=True),
	PyRunner(command='checkPorosity()',iterPeriod=porCheckInterval),
	PyRunner(command='checkPoint()',iterPeriod=checkPointInterval)
]

u='_'

if runningInBatch():
	O.run()
else:
	print('normal run')
	done=True

#This function is required for when the target strain is reached. May be useful if your completion criteria is a ceratian REV size rather than a porosity.
def reachedTargetStrain():
	print('Reached Target Strain')
	O.pause()
	O.wait()

def sphVolume(r):
	return (float(4)/3)*math.pi*r**3

#This function calculates the porosity precisely, assuming there are no spaces where 3 or more spheres occupy the same volume.
def checkPorosity():
	global done,Save
	REVl = O.cell.size
	REVvol = REVl[0]*REVl[1]*REVl[2]
	trueSphVol = 0.0

	#Compute the total volume of all the spheres
	for b in O.bodies:
		sphvol = sphVolume(b.shape.radius)
		trueSphVol += sphvol
	sphCapsVol = 0.0

	#Compute the volume of each of the spherical caps which are overlapping with other spheres
	for i in O.interactions:
		if (i.isReal == True and i.isActive == True and i.geom.penetrationDepth > 0.0):
			R1 = O.bodies[i.id1].shape.radius
			R2 = O.bodies[i.id2].shape.radius
			d = R1 + R2 - i.geom.penetrationDepth
			sphCapsVol += ((math.pi*(R1+R2-d)**2)*(d*d+2*d*R2-3*R2*R2+2*d*R1+6*R2*R1-3*R1*R1))/(12*d)

	#Compute the porosity and exit if criteria is met
	por = (trueSphVol-sphCapsVol)/REVvol
	print('porosity = ', por)
	if(por >= tPorosity):
		print('+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
		print('Target porosity reached')
		print('Finished Simulation')
		print('porosity = ', por)
		print('REVl = ', REVl[0])
		if(Save==1):
			textLoc(REVl, por)
		print('done = ',done)
		done = True
		O.pause()
		O.wait()
	O.engines[6].iterPeriod = max(int(math.ceil(porCheckInterval*(1-(por/tPorosity)**2))),5)
	print('por check interval = ', O.engines[6].iterPeriod)

#Check if the maximum penetration depth is exceeded. If it is larger than Dmax, then the bubble force-displacement interaction equation is being evaluated outside its stated limits and is likely behaving in a non-physical way. See Chan et al. 2010 for details
def checkPen():
	numPen = 0
	for i in O.interactions:
		if (i.isReal == True):
			Ravg = (O.bodies[i.id1].shape.radius+O.bodies[i.id2].shape.radius)/2
			penDepth = i.geom.penetrationDepth
			if (penDepth > i.phys.Dmax):
				numPen += 1
	print('numPen = ', numPen)
	print('fracPen = ', float(numPen)/len(O.interactions))
	return True

def checkUnbalanced():
	ubMaxF = utils.unbalancedForce(True)
	ubAvgF = utils.unbalancedForce(False)
	if ((ubMaxF > alMaxF) | (ubAvgF > alAvgF)):
		print('Unbalanced Force = checkPorosity Fail')
		return False
	else:
		return True

def checkPoint():
	global done
	print('-----------------')
	print('In checkPoint')
	if (checkUnbalanced() & checkPen()):
		print('good run!')
		if(checkPoint.badRun > 0):
			if(checkPoint.badRun > 1):
				O.engines[7].iterPeriod *= 2
			else:
				O.engines[7].iterPeriod = checkPointInterval
			checkPoint.badRun -= 1
#Uncomment this code if you want to try the adaptive time-stepping
		else:
			print('increasing timestep')
			O.engines[5].targetDt *= 1.05
			if(O.engines[5].maxDt < O.engines[5].targetDt):
				O.engines[5].maxDt = O.engines[5].targetDt
		O.saveTmp(str(Row))
		O.pause()
		O.wait()
		O.engines[5].timeStepUpdateInterval = 1
		O.step()
		O.engines[5].timeStepUpdateInterval = O.iter+checkPointInterval+100
		print('O.dt = ',O.dt)
		print('+++++++++++++++++++++++++++++')
		O.run()
	else:
		print('bad run!')
		if (checkPoint.badRun >= maxBadRun):
			sphOut=open(str(Row)+u+'sphin'+'.txt','w')
			sphOut.write('The run failed. Unstable too many times in a row')
			sphOut.close()
			done = True
			O.pause()
			O.wait()
		else:
			checkPoint.badRun += 1
			print('updated badRun = ', checkPoint.badRun)
			print('restarting...')
			import thread
			thread.start_new_thread(restartFunc,())
checkPoint.badRun = 0

def restartFunc():
	tempDt = O.engines[5].targetDt
	tempIter = O.engines[7].iterPeriod
	O.pause()
	O.wait()
	O.loadTmp(str(Row))
	O.engines[7].iterPeriod = tempIter/2
	O.engines[5].timeStepUpdateInterval = 1
	O.engines[5].targetDt = 0.6*tempDt
	O.engines[5].maxDt = O.engines[5].targetDt
	O.step()
	O.engines[5].timeStepUpdateInterval = O.iter+checkPointInterval+100
	O.run()
	print('O.dt = ',O.dt)
	print('+++++++++++++++++++++++++++++')

def writePen(Row):
	global t
	penOut=open(str(Row)+u+'intin'+'.pkl','wb')
	for i in O.interactions:
		penDepth = i.geom.penetrationDepth
		if(penDepth > 0.0):
			Ravg = (O.bodies[i.id1].shape.radius+O.bodies[i.id2].shape.radius)/2
			Fratio = i.phys.fN/(2*math.pi*t*Ravg)
			cPickle.dump(Fratio, penOut)
	penOut.close()

def writeSph(sphOut,crds,r,lim2,lim):
	wrp = [[0,0] for i in range(3)]
	for i in range(3):
		if(crds[i]-r < -lim2[i]):
			wrp[i][1] = 1
		if(crds[i]+r > lim2[i]):
			wrp[i][0] = -1
	
	sCount = 0	
	for i in range(wrp[0][0],wrp[0][1]+1):
		for j in range(wrp[1][0],wrp[1][1]+1):
			for k in range(wrp[2][0],wrp[2][1]+1):
				sphOut.write('\t%g\t%g\t%g\t%g\n'%(crds[0]+i*lim[0],crds[1]+j*lim[1],crds[2]+k*lim[2],r))
				sCount += 1

	return sCount

def textLoc(REVl,por,mask=-1):
	global rMean
	O=Omega()
	writePen(Row)
	R0 = REVl[0]
	R1 = REVl[1]
	R2 = REVl[2]
	try:
		sphOut=open(str(Row) + u + 'sphin.txt','w')
		comOut=open(str(Row) + u + 'comin.txt','w')
	except:
		raise RuntimeError("Problem to write into the file")
	count=0
	extCount=0
	minVec = O.bodies[0].state.pos
	offVec = Vector3(REVl[0]/2,REVl[1]/2,REVl[2]/2)
	bblList = []
	itrList = range(-1,2)
	for b in O.bodies:
		try:
			if(b.state.pos.norm() < minVec.norm()):
				minVec = b.state.pos
		except AttributeError:
			pass
	for b in O.bodies:
		count+=1
		curVec = O.cell.wrap(b.state.pos-minVec)-offVec
		extCount += writeSph(sphOut,curVec,b.shape.radius,offVec,REVl)

	comOut.write('\t%g\n\t%g\n\t%g\n'%(REVl[0],REVl[1],REVl[2]))
	comOut.write('\t%g\n'%(por))
	comOut.write('\t%g\n'%(rMean))
	comOut.write('\t%g\n'%(count))
	comOut.write('\t%g\n'%(extCount))
	comOut.close()
	sphOut.close()
	np.save(str(Row)+u+'sphin',np.array(bblList))

if(c == 0):
	batRet=open('initData.pkl','ab')
	REVl=O.cell.size
	REVvol = REVl[0]*REVl[1]*REVl[2]
	por=sphVol/REVvol
	cPickle.dump(por, batRet)
	cPickle.dump(r, batRet)
	cPickle.dump(numSpheres, batRet)
	batRet.close()

#Don't finish the script unless done=True
while not done:
	utils.waitIfBatch()

