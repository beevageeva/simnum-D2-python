import numpy as np
import sys
from constants import gamma

import initcond_soundwave as initcond

from notifier_params import plotVelFFT

#showErr = True
showErr = False

def getNotifier(notifierType, z, titles, iniValues):
	if(notifierType=="visual"):
		#this can be any other class implementing corresponding methods
		from visual_plot import VisualPlot
		return VisualPlot(z, titles, iniValues)
	elif(notifierType=="file"):
		from save_to_file import SaveToFile 
		return SaveToFile(z, titles, iniValues)

class Model:
	
	def __init__(self):
		from common import getZArray
		from alg import getInitialUcUe
		self.z = getZArray()
		r = initcond.getInitialPresRhoVel(self.z)	
		self.pres = r["pres"]
		self.rho = r["rho"]
		self.vel = r["vel"]
#		print(self.pres.shape)	
#		print("--------initial pres")
#		print(self.pres)
#		print("--------initial rho")
#		print(self.rho)
#		print("--------initial  vel")
#		print(self.vel)
	
		r = getInitialUcUe(self.rho, self.vel, self.pres)
		self.uc = r['uc']
		self.ue = r['ue']
		#flux arrays
		self.fm = None
		self.fc = None
		self.fe = None
		from notifier_params import notifierType
		#self.notifier = getNotifier(notifierType, self.z, ["pres", "rho", "vel"], [self.pres, self.rho, self.vel])
		#self.notifier = getNotifier(notifierType, self.z, ["pres", "rho", "vel", "ue", "uc"], [self.pres, self.rho, self.vel, self.ue, self.uc])
		self.notifier = getNotifier(notifierType, self.z, ["pres", "rho", "vel", "ue", "uc", "ucr"], [self.pres, self.rho, self.vel, self.ue, self.uc, (self.uc[:,:,0]**2 + self.uc[:,:,1]**2)/ (2*self.rho) ])
		if(plotVelFFT):
			self.notifier.addFFTAxis("velFFT0", self.vel[...,0])
			self.notifier.addFFTAxis("velFFT1", self.vel[...,1])
		self.notifier.afterInit()


	def showVars(self):
		#print("x")
		#print(self.z[0])
		#print("y")
		#print(self.z[1])
		print("Complete arrays:")
		np.set_printoptions(threshold='nan')
		print("rho")
		print(self.rho)
		print("pres")
		print(self.pres)
		print("vel")
		print(self.vel)
		print("uc")
		print(self.uc)
		print("ue")
		print(self.ue)
		print("fm")
		print(self.fm)
		print("fc")
		print(self.fc)
		print("fe")
		print(self.fe)
#		print("Proj at 0:")
#		print("rho projx=0")
#		print(self.rho[0,:])
#		print("rho projy=0")
#		print(self.rho[:,0])
#		print("pres projx=0")
#		print(self.pres[0,:])
#		print("pres  projy=0")
#		print(self.pres[:,0])
#		print("vel_x projx=0")
#		print(self.vel[0,:,0])
#		print("vel_x projy=0")
#		print(self.vel[:,0,0])
#		print("vel_y projx=0")
#		print(self.vel[0,:,1])
#		print("vel_y projy=0")
#		print(self.vel[:,0,1])
#		print("uc_x, projx=0")
#		print(self.uc[0,:,0])
#		print("uc_x, projy=0")
#		print(self.uc[:,0,0])
#		print("uc_y, projx=0")
#		print(self.uc[0,:,1])
#		print("uc_y, projy=0")
#		print(self.uc[:,0,1])
#		print("ue projx=0")
#		print(self.ue[0,:])
#		print("ue projy=0")
#		print(self.ue[:,0])
#		print("fm_x, projx=0")
#		print(self.fm[0,:,0])
#		print("fm_x, projy=0")
#		print(self.fm[:,0,0])
#		print("fm_y, projx=0")
#		print(self.fm[0,:,1])
#		print("fm_y, projy=0")
#		print(self.fm[:,0,1])
#		print("fc_xx, projx=0")
#		print(self.fc[0,:,0])
#		print("fc_xx, projy=0")
#		print(self.fc[:,0,0])
#		print("fc_xy, projx=0")
#		print(self.fc[0,:,1])
#		print("fc_xy, projy=0")
#		print(self.fc[:,0,1])
#		print("fc_yy, projx=0")
#		print(self.fc[0,:,2])
#		print("fc_yy, projy=0")
#		print(self.fc[:,0,2])
#		print("fe_x, projx=0")
#		print(self.fe[0,:,0])
#		print("fe_x, projy=0")
#		print(self.fe[:,0,0])
#		print("fe_y, projx=0")
#		print(self.fe[0,:,1])
#		print("fe_y, projy=0")
#		print(self.fe[:,0,1])
		print("END")



	def mainLoop(self, timeEnd):
		import datetime
		from alg import recalculateFluxes, getTimestep, recalculateU, recalculateVelPres,getInitialUcUe
		time = 0.0
		nstep = 0
		ndt = 0
		while(time<timeEnd):
			r = recalculateFluxes(self.rho, self.uc, self.ue, self.vel, self.pres)
			self.fm = r['fm']
			self.fc = r['fc']
			self.fe = r['fe']
			#print("before getTimestep time = %E" % time)
			dt = getTimestep(self.vel, self.pres, self.rho)
			#check if dt is 0 -> because  pres/rho might have become negative see  getTimestep in alg.py
			if dt==0:
				print("STOP dt = 0")
				import time
				time.sleep(5)
				break
			time+=dt
			ndt += dt
			#print("VARS before recalculate at time %4.3f" % time)
			#self.showVars()
			nstep+=1
			#recalculate u at next step - time	
			print("before recalculate at time %4.3f" % time)
			print(datetime.datetime.now())
			result = recalculateU(self.rho, self.uc, self.ue, self.fm, self.fc, self.fe, dt)
			print("after recalculate at time %4.3f" % time)
			print(datetime.datetime.now())
			self.rho = result["rho"]
			self.uc = result["uc"]
			self.ue = result["ue"]
			#recalculate velocity and pression
			r = recalculateVelPres(self.rho, self.uc, self.ue)
			self.vel = r["vel"]	
			self.pres = r["pres"]
			##print("NSTEP %d" % nstep)
			#if time>0.08:
			#	print("VARS after recalculate at time %4.3f" % time)
			#	self.showVars()
			from notifier_params import nstepsPlot, plotAnalitical
			if(nstep % nstepsPlot == 0):
				##print("upd")
				if(plotAnalitical):
					anRes = initcond.getAnRhoPresVel(self.z, time)
					if(showErr):
						err = np.max(np.absolute(np.subtract(self.pres, anRes["pres"])))
						print("time=%E,max abs (self.pres - anPres) = %E" % (time,err))
						err = np.max(np.absolute(np.subtract(self.rho, anRes["rho"])))
						print("max abs (self.rho - anRho) = %E" % err)
						err = np.max(np.absolute(np.subtract(self.vel[:,:,0], anRes["vel"][:,:,0])))
						print("max abs (self.vel0 - anVel0) = %E" % err)
						err = np.max(np.absolute(np.subtract(self.vel[:,:,1], anRes["vel"][:,:,1])))
						print("max abs (self.vel1 - anVel1) = %E" % err)

				#ndt correct value I calculate max velocity only when plotted
				self.notifier.updateValues("rho", [self.rho, anRes["rho"]] if plotAnalitical else self.rho,ndt)
				self.notifier.updateValues("pres", [self.pres, anRes["pres"]] if plotAnalitical else self.pres,ndt)
				self.notifier.updateValues("vel", [self.vel, anRes["vel"]] if plotAnalitical else self.vel,ndt)
				self.notifier.updateValues("ue", self.ue,ndt)
				self.notifier.updateValues("uc", self.uc,ndt)
				self.notifier.updateValues("ucr", (self.uc[:,:,0]**2 + self.uc[:,:,1]**2)/ (2*self.rho) ,ndt)
				if(plotVelFFT):
					self.notifier.updateFFTAxis("velFFT0", self.vel[...,0])
					self.notifier.updateFFTAxis("velFFT1", self.vel[...,1])
				ndt = 0
				self.notifier.afterUpdateValues(time)
		self.notifier.finish()
