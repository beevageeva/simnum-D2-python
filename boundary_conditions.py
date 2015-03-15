import numpy as np

#periodicType = "repeat"  #used for moving plane 
#periodicType = "refl"   #use it with wave packet
periodicType = "diff"  #tried to use it with hankel
		
if periodicType == "repeat":

	def lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		#insert rows	
		array = np.insert(array, 0,  array[n-skip,:], axis = 0)
		array = np.insert(array, n+2,  array[1+skip,:], axis = 0)
		#insert columns	
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n+2, array[:,1+skip], axis = 1)
		return array

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

elif periodicType == "refl" or periodicType == "diff":
		

	def lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		fr = 2 * array[0,:] - array[1,:]
		lr = 2 * array[-1,:] - array[-2,:]
		array = np.insert(array, 0, fr, axis = 0)
		array = np.insert(array, n+2, lr, axis = 0)
		fc = 2 * array[:,0] - array[:,1]
		lc = 2 * array[:,-1] - array[:,-2]
		array = np.insert(array, 0, fc, axis = 1)
		array = np.insert(array, n+2, lc, axis = 1)
		#print("array shape2")
		#print(array.shape)
		#print("end")
		return array


	if periodicType == "refl":
		def lrBoundaryConditionsVel(array, skip=0):
			n = array.shape[0] - 1
			if(skip==0):
				array = np.insert(array, 0,  -array[0,:], axis = 0)
				array = np.insert(array, n+2,  -array[-1,:], axis = 0)
				array = np.insert(array, 0,  -array[:,0], axis = 1)
				array = np.insert(array, n+2,  -array[:,-1], axis = 1)
			elif (skip==1):
				array[0,:] = 0
				array = np.insert(array, 0,  -array[2,:], axis = 0)
				array[-1,:] = 0			
				array = np.insert(array, n+2,  -array[-2,:], axis = 0)
				array[:,0] = 0
				array = np.insert(array, 0,  -array[:,2], axis = 1)
				array[:,-1] = 0			
				array = np.insert(array, n+2,  -array[:,-2], axis = 1)
			return array
	
	elif periodicType == "diff":
		lrBoundaryConditionsVel = lrBoundaryConditionsPresRho
