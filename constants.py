gamma = 5.0/3

#z0 = [3.1, 3.1]
#zf = [7.4, 7.4]

z0 = [0,0]
#zf = [3200, 15000]
zf = [1, 1]


domainType = "centered"   #argType = r expects a domainType = centered
#domainType = "pos00" #use it with d1 and wave packet
#domainType = "unmodified"

if domainType == "centered":
	#centered
	lhalf = 0.5 *  min((zf[1] - z0[1]),(zf[0] - z0[0]) )
	z0[0] = z0[1] = -lhalf
	zf[0] = zf[1] = lhalf


elif domainType == "pos00":
	#0,0
	zf[0] = zf[0] - z0[0]
	zf[1] = zf[1] - z0[1]
	z0[0] = 0
	z0[1] = 0

	

#nint = (512,1024)
nint = (1024,1024)




timeEnd = 1.0 #if not set as program argument it's taken from here


#schemeType = "lf"  # scheme type may be lf(Lax - Fr) or fg (first generation)
schemeType = "fg" 
#schemeType = "fg2" 
#loopType = "python" 
#loopType = "weave" 
loopType = "cython" 
if schemeType == "lf":
	fcfl = 0.99 #use this for lax - fr scheme type
elif schemeType == "fg" or schemeType=="fg2":
	fcfl = 0.97#use this for first generation scheme
	bcStep = "interm"  #in which step to apply boundary conditions
	#bcStep = "final"

