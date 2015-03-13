gamma = 5.0/3

z0_0 = 3.1
zf_0 = 7.4
z0_1 = 3.1
zf_1 = 7.4

#domainType = "centered"   #argType = r expects a domainType = centered
domainType = "pos00" #use it with d1 and wave packet
#domainType = "unmodified"

if domainType == "centered":
	#centered
	lhalf = 0.5 *  min((zf_1 - z0_1),(zf_0 - z0_0) )
	z0_0 = z0_1 = -lhalf
	zf_0 = zf_1 = lhalf


elif domainType == "pos00":
	#0,0
	zf_0 = zf_0 - z0_0
	zf_1 = zf_1 - z0_1
	z0_0 = 0
	z0_1 = 0

	

nint = 64
#nint=128
#nint =  256
#nint = 1024
#nint = 32




timeEnd = 1.0 #if not set as program argument it's taken from here


#schemeType = "lf"  # scheme type may be lf(Lax - Fr) or fg (first generation)
schemeType = "fg" 
#loopType = "python" 
loopType = "weave" 
if schemeType == "lf":
	fcfl = 0.99 #use this for lax - fr scheme type
elif schemeType == "fg":
	fcfl = 0.97#use this for first generation scheme
	bcStep = "interm"  #in which step to apply boundary conditions
	#bcStep = "final"

