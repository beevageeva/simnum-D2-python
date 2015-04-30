import numpy as np
cimport numpy as np

import cython
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_singlestep_u_array_2d(np.ndarray[DTYPE_t, ndim=2] res, np.ndarray[DTYPE_t, ndim=2] u, np.ndarray[DTYPE_t, ndim=3] f, int nint1, int nint2, float dz0, float dz1, float dt):
	cdef int i,j
	for i in range(1, nint1+1):
		for j in range(1, nint2+1):
			#averaging on the first term makes the scheme stable (see Appendix: The lax-Fr scheme)
			res[i-1,j-1] = 0.25 * (u[i,j-1] + u[i,j+1] + u[i+1,j] + u[i-1,j]) - 0.5 * dt * ((f[i+1,j,0] - f[i-1,j,0])/dz0 + (f[i,j+1,1] - f[i,j-1,1])/dz1) 

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_singlestep_u_array_3d(np.ndarray[DTYPE_t, ndim=3] res, np.ndarray[DTYPE_t, ndim=3] u, np.ndarray[DTYPE_t, ndim=3] f, int nint1, int nint2, float dz0, float dz1, float dt):
	cdef int i,j
	for i in range(1, nint1+1):
		for j in range(1, nint2+1):
			res[i-1,j-1,0] = 0.25 * (u[i,j-1,0] + u[i,j+1,0] + u[i+1,j,0] + u[i-1,j,0]) - 0.5 * dt  * ((f[i+1,j,0] - f[i-1,j,0])/dz0 + (f[i,j+1,1] - f[i,j-1,1])/dz1) 
			res[i-1,j-1,1] = 0.25 * (u[i,j-1,1] + u[i,j+1,1] + u[i+1,j,1] + u[i-1,j,1]) - 0.5 * dt * ((f[i+1,j,1] - f[i-1,j,1])/dz0 + (f[i,j+1,2] - f[i,j-1,2])/dz1) 

