from constants import gamma
import numpy as np
from boundary_conditions import lrBoundaryConditionsPresRho, lrBoundaryConditionsVel
lrBoundaryConditions = None


def power2Vec(vfr):
	return vfr[:,:,0] ** 2 + vfr[:,:,1] ** 2


def getInitialUcUe(rho, v , p):
	uc = np.dstack((rho * v[:,:,0], rho * v[:,:,1] )) #cannot directly multiply 2d and 3d array
	ue = p / (gamma - 1.0) + 0.5 * rho * power2Vec(v)
	#print("getInitialUcUe uc_xxPX:")
	#print(" ".join(map(str, uc[0,:,0])))
	#print("getInitialUcUe uePX:")
	#print(" ".join(map(str, ue[0,:])))

	return {'uc': uc, 'ue': ue}

def recalculateVelPres(rho, uc, ue):
	#print("Rho=")
	#print(rho)
	#print("Uc in recalculate Vel Pres ")
	#print(uc)

	v = np.dstack((np.divide(uc[:,:,0], rho ), np.divide(uc[:,:,1], rho )  )) #cannot directly divide 3d to 2d array
	p = (gamma - 1.0) * (ue - power2Vec(uc) / (2.0 * rho) )

	#n = len(p[0])
	#print("max p hor proj at the middle  %e at index %d" % (np.max(p[n/2,:]) , np.argmax(p[n/2,:])) )
	#print("max t1  hor proj at the middle  %e at index %d " %  (np.max(t1[n/2,:]) , np.argmax(t1[n/2,:]))  )
	#print("max v[0]  hor proj at the middle  %e at index %d" %  (np.max(v[n/2,:,0]) , np.argmax(v[n/2,:,0]))  )
#	np.set_printoptions(threshold='nan')
#	print("alg.py: recalculatinf v")
#	print(v)
#	print("alg.py: recalculating pres")
#	print(p)
#	print("uc is ")
#	print(uc)
	#print("vel ")
	#print(v)
	return {'vel': v, 'pres': p}

def recalculateFluxes(rho, uc, ue, v, p):
	fm = uc
	fe1 = ue + p	
	fe = np.dstack((fe1 * v[:,:,0], fe1 * v[:,:,1] )) #cannot directly multiply 2d and 3d array
	fc1 = np.multiply(rho, v[:,:,0] ** 2) + p
	fc2 = np.multiply(rho, np.multiply(v[:,:,0], v[:,:,1])) #Fc_xz = Fc_zx store only ince the value
	fc3 = np.multiply(rho, v[:,:,1] ** 2) + p
	#print("recalculateFluxes fc3")
	#print(fc3)	
	#print("recalculateFluxes v_z")
	#print(v[:,:,1])	
	#print("recalculateFluxes v_z**2")
	#print(v[:,:,1] ** 2)
	#!!! dstack won't work as expected for 3 ? it's the same	
	fc = np.dstack((fc1,fc2,fc3))
	#fc =  np.zeros(fc1.shape + (3,), dtype='float')
	#fc[:,:, 0] = fc1
	#fc[:,:, 1] = fc2
	#fc[:,:, 2] = fc3
	#print("recalculateFluxes pres")
	#print(p)
	#print("recalculateFluxes fc")
	#print(fc)
	#print("recalculateFluxes fc third column")
	#print(fc[:,:,2])
	return {'fm': fm, 'fc': fc, 'fe': fe}

def getTimestep(v, p, rho):
	
	#print("getTimestep presPX")
	#print(" ".join(map(str, p[0,:])))
	#print("getTimestep rhoPX")
	#print(" ".join(map(str, rho[0,:])))
	#print("getTimestep vel_xPX")
	#print(" ".join(map(str, v[0,:,0])))
	from constants import fcfl
	from common import getDz0, getDz1
	dz0 = getDz0()
	dz1 = getDz1()
	t1 =  np.divide(p, rho)
	if(np.any(t1<0)):
		print("there is p/rho < 0 in getTimestep")
		#print("p<0")
		#print(p[p<0])
		#print("rho<0")
		#print(rho[rho<0])
		print(np.where(t1 < 0)[0])

		return 0	
	cs = np.sqrt(gamma *  t1)
	#the extension of 1D is unstable because on Neumann critrion
	#smax0 = np.max(np.concatenate([np.absolute(v[:,:,0] + cs), np.absolute(v[:,:,0] - cs)]))
	#smax1 = np.max(np.concatenate([np.absolute(v[:,:,1] + cs), np.absolute(v[:,:,1] - cs)]))
	#dt = 0.5 * min(float(dz0 * fcfl ) / smax0, float(dz1 * fcfl ) / smax1 )
	#see the following link
	#TODO	
	#http://homepage.univie.ac.at/franz.vesely/cp_tut/nol2h/new/c5pd_s1ih.html 
	#http://www.aei.mpg.de/~rezzolla/lnotes/Evolution_Pdes/evolution_pdes_lnotes.pdf pag 59
	smax = np.max(np.concatenate([np.sqrt( (v[:,:,0] + cs) ** 2 + (v[:,:,1] + cs)**2 ), np.sqrt( (v[:,:,0] - cs) ** 2 + (v[:,:,1] - cs)**2 )]))
	dt = min(dz0, dz1) * fcfl  / (smax *  2 ** (0.5))
#	#print("getTimestep %E" % dt)

	return dt



from constants import schemeType, loopType

if loopType == "cython":
	import pyximport
	pyximport.install(setup_args={'include_dirs':[np.get_include()]})


	
if schemeType == "fg":

	def calcIntermUArray(u, f, dt):
		#print("calcIntermU")
		from common import getDz0, getDz1
		from constants import nint
		#no more lambda dz1 may be different from dz0
		dz0 = getDz0()
		dz1 = getDz1()
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((nint+1, nint+1))
			#res = np.array(nint+1, nint+1) #TODO do not initialize how?

			if loopType == "python":
				from python_alg2_fg import calc_interm_u_array_2d
			elif loopType == "weave":
				from weave_alg2_fg import calc_interm_u_array_2d
			elif loopType == "cython":
				from cython_alg2_fg import calc_interm_u_array_2d
			calc_interm_u_array_2d(res, u, f, nint, dz0, dz1, dt) 


		else:
			res = np.zeros((nint+1, nint+1, 2))
			if loopType == "python":
				from python_alg2_fg import calc_interm_u_array_3d
			elif loopType == "weave":
				from weave_alg2_fg import calc_interm_u_array_3d
			elif loopType == "cython":
				from cython_alg2_fg import calc_interm_u_array_3d
			calc_interm_u_array_3d(res, u, f, nint, dz0, dz1, dt) 
	
		return res

	
	def calcFinalUArray(u, intermF, dt, skip=0):
		#print("calcFinalU")
		from common import getDz0, getDz1
		from constants import nint
		dz0 = getDz0()
		dz1 = getDz1()
		n = intermF.shape[0] - 1
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((n, n))
			#res = np.array(nint+1, nint+1) #TODO do not initialize how?

			if loopType == "python":
				from python_alg2_fg import calc_final_u_array_2d
			elif loopType == "weave":
				from weave_alg2_fg import calc_final_u_array_2d
			elif loopType == "cython":
				from cython_alg2_fg import calc_final_u_array_2d
			calc_final_u_array_2d(res, u, intermF, n, dz0, dz1, dt, skip) 

		else:
			res = np.zeros((n, n, 2))

			if loopType == "python":
				from python_alg2_fg import calc_final_u_array_3d
			elif loopType == "weave":
				from weave_alg2_fg import calc_final_u_array_3d
			elif loopType == "cython":
				from cython_alg2_fg import calc_final_u_array_3d
			calc_final_u_array_3d(res, u, intermF, n, dz0, dz1, dt, skip) 

		#no more boundary conditions because intermediate array alreday has nint + 3 points
		return res

	from constants import bcStep
	if(bcStep == "interm"):
		def calcIntermU(u, f, dt):
			res = calcIntermUArray(u, f, dt)
			#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 3 points see array limits
			#print("calcIntermStep before bc")
			#print(res)	
			res = lrBoundaryConditions(res, 1)
			#print("calcIntermStep after bc")
			#print(res)	
			return res
	
		#no more boundary conditions because intermediate array alreday has nint + 3 points
		calcFinalU = calcFinalUArray	
		
	elif(bcStep == "final"):
		#boundary conditions in final step
			#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 1 points see array limits
		calcIntermU = calcIntermUArray
		def calcFinalU(u, f, dt):
			res = calcFinalUArray(u, f, dt,1)
			#print("final bc array before bc")
			#print(res)
			res = lrBoundaryConditions(res)
			#print("final bc array after bc")
			#print(res)
			return res
	
	def recalculateU(rho, uc, ue, fm, fc ,fe, dt):
		global lrBoundaryConditions
		lrBoundaryConditions = lrBoundaryConditionsPresRho
		#print("recalculateU ")
		intermRho = calcIntermU(rho, fm , dt)	
		intermUe = calcIntermU(ue, fe , dt)
		lrBoundaryConditions = lrBoundaryConditionsVel
		intermUc = calcIntermU(uc, fc , dt)	
		intermVelPres = recalculateVelPres(intermRho, intermUc, intermUe)
		intermVel = intermVelPres["vel"]
		intermPres = intermVelPres["pres"]
		intermFluxes = recalculateFluxes(intermRho, intermUc, intermUe, intermVel, intermPres)
		lrBoundaryConditions = lrBoundaryConditionsPresRho
		finalRho = calcFinalU(rho, intermFluxes["fm"], dt)	
		finalUe = calcFinalU(ue, intermFluxes["fe"], dt)
		lrBoundaryConditions = lrBoundaryConditionsVel
		finalUc = calcFinalU(uc, intermFluxes["fc"], dt)	
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}


elif schemeType == "fg2":
	def calcIntermUArray(u, f, dt):
		#print("calcIntermU")
		from common import getDz0, getDz1
		from constants import nint
		#no more lambda dz1 may be different from dz0
		dz0 = getDz0()
		dz1 = getDz1()
		if(u.ndim == 2): #case of um and ue that are not vectors
			resx = np.zeros((nint+1, nint))
			resy = np.zeros((nint, nint+1))
			#res = np.array(nint+1, nint+1) #TODO do not initialize how?
			if loopType == "python":
				from python_alg2_fg2 import calc_interm_u_array_2d
			elif loopType == "weave":
				from weave_alg2_fg2 import calc_interm_u_array_2d
			elif loopType == "cython":
				from cython_alg2_fg2 import calc_interm_u_array_2d
				
			calc_interm_u_array_2d(resx, resy, u,f,nint, dz0, dz1, dt) 

		else:
			resx = np.zeros((nint+1, nint, 2))
			resy = np.zeros((nint, nint+1, 2))
			#res = np.array(nint+1, nint+1, 2)
			if loopType == "python":
				from python_alg2_fg2 import calc_interm_u_array_3d
			elif loopType == "weave":
				from weave_alg2_fg2 import calc_interm_u_array_3d
			elif loopType == "cython":
				from cython_alg2_fg2 import calc_interm_u_array_3d
			calc_interm_u_array_3d(resx, resy, u, f, nint, dz0, dz1, dt) 
	
		return resx, resy

	
	def calcFinalUArray(u, intermF0, intermF1, dt, skip=0):
		#print("calcFinalU")
		from common import getDz0, getDz1
		from constants import nint
		dz0 = getDz0()
		dz1 = getDz1()
		n = intermF0.shape[1]  #use as shape index the dimension d which is = other dim - 1
		print("final fg2 array n = %d" % n) 
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((n, n))
			#res = np.array(nint+1, nint+1) #TODO do not initialize how?
			if loopType == "python":
				from python_alg2_fg2 import calc_final_u_array_2d
			elif loopType == "weave":
				from weave_alg2_fg2 import calc_final_u_array_2d
			elif loopType == "cython":
				from cython_alg2_fg2 import calc_final_u_array_2d
				calc_final_u_array_2d(res, u, intermF0, intermF1, n, dz0, dz1, dt, skip) 
		else:
			res = np.zeros((n, n, 2))
			#res = np.array(nint+1, nint+1, 2)
			if loopType == "python":
				from python_alg2_fg2 import calc_final_u_array_3d
			elif loopType == "weave":
				from weave_alg2_fg2 import calc_final_u_array_3d
			elif loopType == "cython":
				from cython_alg2_fg2 import calc_final_u_array_3d
			calc_final_u_array_3d(res, u, intermF0, intermF1, n, dz0, dz1, dt, skip) 
		#no more boundary conditions because intermediate array alreday has nint + 3 points
		return np.array(res)


	from constants import bcStep
	if(bcStep == "interm"):
		def calcIntermU(u, f, dt):
			resx, resy = calcIntermUArray(u, f, dt)
			#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 3 points see array limits
			#print("calcIntermStep before bc")
			#print(res)	
			resx = lrBoundaryConditions(resx, 2)
			resy = lrBoundaryConditions(resy, 2)
			#print("calcIntermStep after bc")
			#print(res)	
			return resx, resy
	
		#no more boundary conditions because intermediate array alreday has nint + 3 points
		calcFinalU = calcFinalUArray	
		
	elif(bcStep == "final"):
		#boundary conditions in final step
			#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 1 points see array limits
		calcIntermU = calcIntermUArray
		def calcFinalU(u, fx, fy, dt):
			res = calcFinalUArray(u, fx, fy, dt,1)
			#print("final bc array before bc")
			#print(res)
			res = lrBoundaryConditions(res)
			#print("final bc array after bc")
			#print(res)
			return res
	
	def recalculateU(rho, uc, ue, fm, fc ,fe, dt):
		global lrBoundaryConditions
		lrBoundaryConditions = lrBoundaryConditionsPresRho
		#print("recalculateU ")
		print("------------before interm u 2d----------")
		intermRho0, intermRho1 = calcIntermU(rho, fm , dt)	
		#print("------------intermrho0----------")
		#print(intermRho0)
		#print(intermRho1)
		#print("------------irho0end----------")
		intermUe0, intermUe1 = calcIntermU(ue, fe , dt)
		lrBoundaryConditions = lrBoundaryConditionsVel
		print("------------before interm u 3d----------")
		intermUc0, intermUc1 = calcIntermU(uc, fc , dt)	
#		print("------------INTERM SHAPES0----------")
#		print(intermRho0.shape)
#		print(intermUe0.shape)
#		print(intermUc0.shape)
#		print("------------INTERM SHAPES0-END---------")
#		print("------------INTERM SHAPES1----------")
#		print(intermRho1.shape)
#		print(intermUe1.shape)
#		print(intermUc1.shape)
#		print("------------INTERM SHAPES1-END---------")
		intermVelPres = recalculateVelPres(intermRho0, intermUc0, intermUe0)
		intermVel0 = intermVelPres["vel"]
		intermPres0 = intermVelPres["pres"]
		intermVelPres = recalculateVelPres(intermRho1, intermUc1, intermUe1)
		intermVel1 = intermVelPres["vel"]
		intermPres1 = intermVelPres["pres"]
		intermFluxes0 = recalculateFluxes(intermRho0, intermUc0, intermUe0, intermVel0, intermPres0)
		intermFluxes1 = recalculateFluxes(intermRho1, intermUc1, intermUe1, intermVel1, intermPres1)
		lrBoundaryConditions = lrBoundaryConditionsPresRho
		print("------------before final u 2d----------")
		finalRho = calcFinalU(rho, intermFluxes0["fm"], intermFluxes1["fm"], dt)	
		finalUe = calcFinalU(ue, intermFluxes0["fe"], intermFluxes1["fe"], dt)
		lrBoundaryConditions = lrBoundaryConditionsVel
		print("------------before final u 3d----------")
		finalUc = calcFinalU(uc, intermFluxes0["fc"], intermFluxes1["fc"], dt)	
#		print("------------FINAL SHAPES----------")
#		print(finalRho.shape)
#		print(finalUc.shape)
#		print(finalUe.shape)
#		print("------------FINAL SHAPES-END---------")

		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}

		
	

elif schemeType == "lf":

	def calcSingleStepU(u,f, dt):
		#print("calcSingleStepU f first comp")
		#print(f[:,:,0])
		#print("calcSingleStepU f second comp")
		#print(f[:,:,1])
		from common import getDz0, getDz1
		from constants import nint
		dz0 = getDz0()
		dz1 = getDz1()
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((nint, nint))
			if loopType == "python":
				from python_alg2_lf import calc_singlestep_u_array_2d
			elif loopType == "weave":
				from weave_alg2_lf import calc_singlestep_u_array_2d
			elif loopType == "cython":
				from cython_alg2_lf import calc_singlestep_u_array_2d
			calc_singlestep_u_array_2d(res, u,f,nint, dz0,dz1, dt) 
		else:
			res = np.zeros((nint, nint, 2))
			if loopType == "python":
				from python_alg2_lf import calc_singlestep_u_array_3d
			elif loopType == "weave":
				from weave_alg2_lf import calc_singlestep_u_array_3d
			elif loopType == "cython":
				from cython_alg2_lf import calc_singlestep_u_array_3d
			calc_singlestep_u_array_3d(res, u,f,nint, dz0,dz1, dt) 
		#print("calcSingleStep before BC PX")
		#print(res[0,...])
		res = lrBoundaryConditions(res)
		#print("calcSingleStep after BC PX")
		#print(res[0,...])
		#print("res AFTER BC")
		#print(res)
		return res
	
	def recalculateU(rho, uc, ue, fm , fc, fe, dt):
		global lrBoundaryConditions
		lrBoundaryConditions = lrBoundaryConditionsPresRho
		#print("recalculateRho")
		finalRho = calcSingleStepU(rho, fm, dt)
		#print("recalculateUe")
		finalUe = calcSingleStepU(ue, fe, dt)
		lrBoundaryConditions = lrBoundaryConditionsVel
		#print("uc  (before recalculateU) = ")
		#print(uc)	
		#print("fc (recalculateU) = ")
		#print(fc)	
		#print("fc third column (recalculateU) = ")
		#print(fc[:,:,2])	
		#print("recalculateUc")
		finalUc = calcSingleStepU(uc, fc, dt)	
		#print("uc  (after recalculateU) = ")
		#print(finalUc)	
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}

else:
	puts("Scheme type not implemented %s" % schemeType)



