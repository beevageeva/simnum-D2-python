import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter


matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys, os
from scipy.fftpack import fft,fftfreq#forFourierTransform

from notifier_params import plots, plotAnalitical

saveImages = True
#saveImages = False



plot2dModule = True
#plot2dModule = False

#ylim = {"pres":{ "maxY": 1.0005, "minY": 0.9995} , "vel" : { "maxY": 0.00035, "minY": -0.00035}, "rho":{ "maxY": 1.0004, "minY": 0.9996}} 
ylim = {"pres":{ "maxY": 1.0003, "minY": 0.9997} , "vel0" : { "maxY": 0.00025, "minY": -0.00025}, "vel1" : { "maxY": 0.00025, "minY": -0.00025}, "rho":{ "maxY": 1.0003, "minY": 0.9997}, "vel" : { "maxY": 0.0005, "minY": -0.0001}} 
#ylim = {"pres":{ "maxY": 1.0003, "minY": 0.9997} , "vel" : { "maxY": 0.00025, "minY": -0.00025}, "rho":{ "maxY": 0.5002, "minY": 0.4998}} 
#ylim = {"pres":{ "maxY": 1.0003, "minY": 0.9997} , "vel" : { "maxY": 0.00025, "minY": -0.00025}, "rho":{ "maxY": 2.0002, "minY": 1.9998}} 
#xlim = {"minX" : 0, "maxX" : 4.3}
#ylim = None
xlim = None

	
#because methods are different in python3 and python2
def testKeyInDict(key, dictionary):
	import sys	
	if (sys.version_info[0]==2):
		return dictionary.has_key(key)
	else:
		return key in dictionary	





class VisualPlot:

	#TODO this does not use self
	def addAxisColor(self, arrayToAppendAxes, title, vals, n , i, subplotNumber, colspan):
		if(vals.ndim == 2):
			values = vals
		else:
			values = vals[..., subplotNumber]
		if(colspan):
			ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=2)
		else:
			ax = plt.subplot2grid((n,2), (i,subplotNumber))
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_title("%s%d"% (title, subplotNumber))
		ax.grid(True)
		#np.set_printoptions(threshold='nan')
		#print("addAxisColor: values")
		#print(values)
		ax.imshow(values)
		ax.relim()
		ax.autoscale_view(True,True,True)
		arrayToAppendAxes.append(ax)

	def addProjAxes(self, arrayToAppendAxes, title, vals, n , i, subplotNumber, colspan=False):
		if(vals.ndim == 2):
			newtitle = title
		else:
			newtitle = "%s%d" % (title, subplotNumber)
		if hasattr(self, "dim0ProjIndex"):
			plt.figure(2)
			if(colspan):
				ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=2)
			else:
				ax = plt.subplot2grid((n,2), (i,subplotNumber))
			if(vals.ndim == 2):
				values = vals[self.dim0ProjIndex, :]
			else:
				values = vals[self.dim0ProjIndex, :, subplotNumber]
			markMaxValue = None
			if(plots["dim0"][1]):
				markMaxIndex = np.argmax(values)
				markMaxValue = self.z[0][0][markMaxIndex]
				#TODO if it works replace this by newtitle
				self.maxPoints["dim0"]["%s%d" % (title, subplotNumber)] = markMaxValue
			self.addAxisProj(ax, newtitle, values, markMaxValue)
			arrayToAppendAxes.append(ax)
		if hasattr(self, "dim1ProjIndex"):
			if(vals.ndim == 2):
				values = vals[:,self.dim1ProjIndex]
			else:
				values = vals[:,self.dim1ProjIndex, subplotNumber]
			plt.figure(3)
			if(colspan):
				ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=2)
			else:
				ax = plt.subplot2grid((n,2), (i,subplotNumber))
			markMaxValue = None
			if(plots["dim1"][1]):
				markMaxIndex = np.argmax(values)
				markMaxValue = self.z[0][0][markMaxIndex]
				self.maxPoints["dim1"]["%s%d" % (title, subplotNumber)] = markMaxValue
			self.addAxisProj(ax, newtitle, values, markMaxValue)
			arrayToAppendAxes.append(ax)

		if hasattr(self, "projMask"):
			if(vals.ndim == 2):
				values = vals[self.projMask]
			else:
				values = vals[self.projMask, subplotNumber]
			plt.figure(6)
			if(colspan):
				ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=2)
			else:
				ax = plt.subplot2grid((n,2), (i,subplotNumber))
			markMaxValue = None
			linez = np.dstack((self.z[0][self.projMask], self.z[1][self.projMask]))
			linez = linez[0]	
			newz = np.sqrt(  (linez[:,0] - linez[0][0])**2 + (linez[:,1] - linez[0][1])**2  )
#			print("addProjAxis linez")
#			print(linez)
#			print("z[0][0]")
#			print(self.z[0][0])
#			print(self.z[0][0].shape)
#			print("newz")
#			print(newz)
#			print(newz.shape)
			if(plots["line"][1]):
				markMaxIndex = np.argmax(values)
				markMaxValue = newz[markMaxIndex]
				self.maxPoints["line"]["%s%d" % (title, subplotNumber)] = markMaxValue
			self.addAxisProj(ax, newtitle, values, markMaxValue, newz)
			arrayToAppendAxes.append(ax)

		if testKeyInDict("color", plots):
			plt.figure(4)
			self.addAxisColor(arrayToAppendAxes, title, vals, n , i, subplotNumber,colspan)
			plt.draw()
			if plotAnalitical:
				plt.figure(5)
				self.addAxisColor(arrayToAppendAxes, ("%s-an" % title), vals, n , i, subplotNumber,colspan)



	
	def addAxisProj(self, ax, title, vals, markMaxValue = None, newz = None):
		if newz is None:
			newz = self.z[0][0]
		ax.set_xlabel("z")
		ax.set_ylabel(title)
		ax.grid(True)
		ax.plot(newz, vals)
		if(not markMaxValue is None):
			ax.vlines(markMaxValue, np.min(vals[0]) if plotAnalitical else np.min(vals), np.max(vals[0]) if plotAnalitical else np.max(vals), color='b', label="max")
			#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		if(ylim and testKeyInDict(title, ylim)):
			ax.set_ylim(ylim[title]["minY"],ylim[title]["maxY"])
			ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
			if xlim:
				ax.set_xlim(xlim["minX"],xlim["maxX"])
			else:
				ax.autoscale_view(True,True,False)
			
		else:
			ax.relim()
			ax.autoscale_view(True,True,True)


	#titles will be an array of  of 2 elem array : the first the title and the second if we will plot 
	def __init__(self, z, titles, iniValues):
		if saveImages:
			from common import createFolder	
			self.dirname = createFolder("outImages")
		#if for example velocity has 2 components I have to have 2 subplots
		self.z = z
		self.maxPoints = {}
		self.figures = []
		if(plots["3d"]):
			fig = plt.figure(1)
			fig.suptitle("3D plot")
			self.figures.append(fig)
		if plots:
			from common import getZIndex0, getZIndex1
			if testKeyInDict("dim0", plots):
				if plots["dim0"][1]:
					self.maxPoints["dim0"] = {}
				fig = plt.figure(2)
				self.dim0ProjIndex = getZIndex0(plots["dim0"][0])
				fig.suptitle("Dim0 z1=%4.3f" % plots["dim0"][0])
				plt.get_current_fig_manager().window.wm_geometry("1000x900+50+50")
				fig.subplots_adjust(right=0.8)
				self.figures.append(fig)
			if testKeyInDict("dim1", plots):
				if plots["dim1"][1]:
					self.maxPoints["dim1"] = {}
				fig = plt.figure(3)
				fig.suptitle("Dim1 z0=%4.3f" % plots["dim1"][0])
				self.dim1ProjIndex = getZIndex1(plots["dim1"][0])
				plt.get_current_fig_manager().window.wm_geometry("1000x900+50+50")
				fig.subplots_adjust(right=0.8)
				self.figures.append(fig)
			if testKeyInDict("line", plots):
				if plots["line"][1]:
					self.maxPoints["line"] = {}
				fig = plt.figure(6)
				fig.suptitle("Line %s" % str(plots["line"][0]))
				n = len(z[0])
				y,x=np.ogrid[-n / 2: n/2 , -n / 2: n/2 ]
				self.projMask = plots["line"][0](x,y)
#				print("************visual_plot PROJMASK")
#				print(self.projMask)
#				print("shape")
#				print(self.projMask.shape)
#				print("************visual_plot z[0]")
#				print(z[0])
#				print("shape")
#				print(z[0].shape)

				plt.get_current_fig_manager().window.wm_geometry("1000x900+50+50")
				fig.subplots_adjust(right=0.8)
				self.figures.append(fig)
			if testKeyInDict("color", plots):
				fig = plt.figure(4)
				fig.suptitle("color")
				self.figures.append(fig)
				if(plotAnalitical):
					fig = plt.figure(5)
					fig.suptitle("color an")
					self.figures.append(fig)
					
		self.axes = {}
		n = len(titles)
		if(plot2dModule):	
			for i in range(0, len(titles)):
				if(iniValues[i].ndim == 3):	
					n+=1
		ii=0
		for i in range(0, len(titles)):
			vals = iniValues[i]
			title = titles[i]
			if(vals.ndim == 2):	
				self.axes[title] = []
				if(plots["3d"]):
					plt.figure(1)
					ax = plt.subplot2grid((n,2), (ii,0), colspan=2, projection='3d')
					self.addAxis(ax, title, vals)
					self.axes[title].append(ax)
				self.addProjAxes(self.axes[title], title, vals, n, ii, 0,  True)
	

			elif(vals.ndim == 3):	
				self.axes[title] = [[]]
				if(plots["3d"]):
					plt.figure(1)
					ax = plt.subplot2grid((n,2), (ii,0), projection='3d')

					self.addAxis(ax, ("%s dim 0" % title), vals[:,:,0])
					self.axes[title][0].append(ax)
				self.addProjAxes(self.axes[title][0], title, vals, n, ii, 0)
				self.axes[title].append([])
				if(plots["3d"]):
					plt.figure(1)
					ax2 = plt.subplot2grid((n,2), (ii,1), projection='3d')
					self.addAxis(ax2, ("%s dim 1" % title), vals[:,:,1])
					self.axes[title][1].append(ax2)
				self.addProjAxes(self.axes[title][1], title, vals, n, ii, 1)
				#module
				if plot2dModule:
					ii+=1
					self.axes[title].append([])
					modVals = np.sqrt(vals[:,:,0] ** 2 + vals[:,:,1] ** 2)
					if(plots["3d"]):
						plt.figure(1)
						ax3 = plt.subplot2grid((n,2), (ii,0),colspan=2, projection='3d')
						self.addAxis(ax3, ("%s mod" % title), modVals)
						self.axes[title][2].append(ax3)
					self.addProjAxes(self.axes[title][2], title, modVals, n, ii, 0, True)


			else:
				#print("Dim invalid %d" % vals.ndim)
				sys.exit(0)
			ii+=1
		if(plots["3d"]):
			plt.figure(1)
			wm = plt.get_current_fig_manager()
			wm.window.wm_geometry("800x900+50+50")
		plt.draw()
		plt.show(block=False)

	def afterInit(self):
		#import time
		#time.sleep(20)
		#save initial figures to files
		if saveImages:
			numFig = 0
			for fig in self.figures:
				os.mkdir(os.path.join(self.dirname, "Fig%d" % numFig))
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img000000.png"))
				numFig +=1

		
	def addAxis(self, ax, title, vals):
		ax.set_xlabel("z")
		ax.set_ylabel("z1")
		ax.set_zlabel(title)
		ax.grid(True)
		#ax.plot_surface(self.z[0], self.z[1], vals, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		#ax.plot_surface(self.z[0], self.z[1], vals)
		ax.plot_wireframe(self.z[0], self.z[1], vals)
		#ax.view_init(0, 90)
		#ax.view_init(45, 45)
		ax.relim()
		ax.autoscale_view(True,True,True)
	

	def updateAxis(self, ax, title, vals):	
		ax.cla()
		ax.set_xlabel("z0")
		ax.set_ylabel("z1")
		ax.set_zlabel(title)
		#no plot on 3d axis for analitical 
		if(plotAnalitical):
			vals = vals[0]
		#ax.plot_surface(self.z[0], self.z[1], vals, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		#ax.plot_surface(self.z[0], self.z[1], vals)
		#import datetime
		#print("before wireframe")
		#print(datetime.datetime.now())
		ax.plot_wireframe(self.z[0], self.z[1], vals)
		#print("after wireframe")
		#print(datetime.datetime.now())
		ax.grid(True)
		ax.relim()
		ax.autoscale_view(True,True,True)

	def updateAxisProj(self, ax, title, vals, markMaxValue = None, maxLegend = None, newz = None):
		if newz is None:
			newz =  self.z[0][0]	
		ax.cla()
		ax.set_xlabel("z")
		ax.set_ylabel(title)
		if(plotAnalitical):
			ax.plot(newz,  vals[0], label="num")
			ax.plot(newz,  vals[1], color="r", label="an")
		else:
#			print("updateAxisProj: title %s " % title)	
#			print("newz")
#			print(newz)
#			print("newz shape")
#			print(newz.shape)
#			print("vals")
#			print(vals)
#			print("vals shape")
#			print(vals.shape)
			ax.plot(newz,  vals)
		if(markMaxValue):
			ax.vlines(markMaxValue, np.min(vals[0]) if plotAnalitical else np.min(vals), np.max(vals[0]) if plotAnalitical else np.max(vals), color='b', label=maxLegend)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		ax.grid(True)
		if(ylim and testKeyInDict(title, ylim)):
			ax.set_ylim(ylim[title]["minY"],ylim[title]["maxY"])
			ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
			if xlim:
				ax.set_xlim(xlim["minX"],xlim["maxX"])
			else:
				ax.autoscale_view(True,True,False)
		else:
			ax.relim()
			ax.autoscale_view(True,True,True)


	def updateAxisColor(self, ax, title, vals):	
		ax.cla()
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_title(title)
		ax.imshow(vals)

	#TODO make a function for every projection
	def updateProjAxis(self, axesArray, title, vals, index, dt):
		from common import getSpeedPeriodic0, getSpeedPeriodic1
		vdim = vals[0].ndim if plotAnalitical else vals.ndim
		if vdim == 2:
			newtitle = title
		else:
			newtitle = "%s%d" % (title, index)
		if plots["3d"]:
			ni=1
		else:
			ni=0
		if hasattr(self, "dim0ProjIndex"):
			if(vdim == 2):
				values = [vals[0][self.dim0ProjIndex, :] , vals[1][self.dim0ProjIndex, :]] if plotAnalitical else vals[self.dim0ProjIndex, :]
			else:
				values = [vals[0][self.dim0ProjIndex, :, index],vals[1][self.dim0ProjIndex, :, index] ]	if plotAnalitical else vals[self.dim0ProjIndex, :, index]	
			markMaxValue = None
			markMaxTitle = None
			if(plots["dim0"][1]):
				markMaxIndex = np.argmax(values[0] if plotAnalitical else values)
				markMaxValue = self.z[0][0][markMaxIndex]
				maxSpeed  = getSpeedPeriodic0(markMaxValue, self.maxPoints["dim0"]["%s%d" % (title, index)], dt)
				markMaxTitle = "ms= %4.3f" % maxSpeed
				self.maxPoints["dim0"]["%s%d" % (title, index)] = markMaxValue
			self.updateAxisProj(axesArray[ni], newtitle, values, markMaxValue, markMaxTitle)
			ni+=1
		if hasattr(self, "dim1ProjIndex"):
			if(vdim == 2):
				values =  [vals[0][:,self.dim1ProjIndex], vals[1][:,self.dim1ProjIndex]] if plotAnalitical else vals[:,self.dim1ProjIndex]
			else:
				values = [vals[0][:,self.dim1ProjIndex, index], vals[1][:,self.dim1ProjIndex, index]]  if plotAnalitical else vals[:,self.dim1ProjIndex, index]
			markMaxValue = None
			markMaxTitle = None
			if(plots["dim1"][1]):
				markMaxIndex = np.argmax(values[0]) if plotAnalitical else np.argmax(values)
				markMaxValue = self.z[0][0][markMaxIndex]
				maxSpeed = getSpeedPeriodic1(markMaxValue, self.maxPoints["dim1"]["%s%d" % (title, index)], dt)	
				markMaxTitle = "ms= %4.3f" % maxSpeed
				self.maxPoints["dim1"]["%s%d" % (title, index)] = markMaxValue
			self.updateAxisProj(axesArray[ni], newtitle, values, markMaxValue, markMaxTitle)
			ni+=1

		if hasattr(self, "projMask"):
			if(vdim == 2):
				values =  [vals[0][self.projMask], vals[1][self.projMask]] if plotAnalitical else vals[self.projMask]
#				if not plotAnalitical:	
#					print("MASK UPDATE PROJ AXIS")
#					print(np.max(vals))
#					print("max on values proj")
#					print(np.max(values))
			else:
				values = [vals[0][:,:,index][self.projMask], vals[1][:,:,index][self.projMask]]  if plotAnalitical else vals[:,:,index][self.projMask]
			markMaxValue = None
			markMaxTitle = None
			linez = np.dstack((self.z[0][self.projMask], self.z[1][self.projMask]))
			linez = linez[0]	
			newz = np.sqrt(  (linez[:,0] - linez[0][0])**2 + (linez[:,1] - linez[0][1])**2  )
			
			if(plots["line"][1]):
				markMaxIndex = np.argmax(values[0]) if plotAnalitical else np.argmax(values)
				markMaxValue = newz[markMaxIndex]
				#TODO speed periodic for line
				#maxSpeed = getSpeedPeriodic1(markMaxValue, self.maxPoints["line"]["%s%d" % (title, index)], dt)	
				#markMaxTitle = "ms= %4.3f" % maxSpeed
				self.maxPoints["line"]["%s%d" % (title, index)] = markMaxValue
			self.updateAxisProj(axesArray[ni], newtitle, values, markMaxValue, markMaxTitle, newz)
			ni+=1

		if testKeyInDict("color", plots):
			if(vdim == 2):
				values = vals[0] if plotAnalitical else vals
			else:
				values = vals[0][:,:,index] if plotAnalitical else vals[:,:, index]
			self.updateAxisColor(axesArray[ni], newtitle, values)
			if(plotAnalitical):
				values = vals[1] if vdim == 2 else vals[1][:,:,index]
				self.updateAxisColor(axesArray[ni+1], newtitle + "-an", values)


	def updateValues(self, title, vals,dt):
		#print("vals.shape")
		#print(vals.shape)

		#print("updateValues %s MAX POINTS:" % title)
		#print(self.maxPoints)
		#print("MAX POINTS END")
		#print("updateValues %s vals.ndim = %d, vals.shape" % (title, vals.ndim))
		#print(vals.shape)
		
		ax = self.axes[title]
		
		vdim = vals[0].ndim if plotAnalitical else vals.ndim
		if(vdim == 2):
			if(plots["3d"]):
				self.updateAxis(ax[0], title, [vals[0], vals[1]] if plotAnalitical else vals)
			self.updateProjAxis(ax, title, vals,0,dt)
		else:
			if plot2dModule:
				if plotAnalitical:
					modVals = [np.sqrt(vals[0][:,:,0] ** 2 + vals[0][:,:,1] ** 2),np.sqrt(vals[1][:,:,0] ** 2 + vals[1][:,:,1] ** 2)]
				else:
					modVals = np.sqrt(vals[:,:,0] ** 2 + vals[:,:,1] ** 2)

			if(plots["3d"]):
				self.updateAxis(ax[0][0],  title, [vals[0][...,0], vals[1][...,0]] if plotAnalitical else vals[...,0])
				self.updateAxis(ax[1][0], title, [vals[0][...,1], vals[1][...,1]] if plotAnalitical else  vals[...,1])
				if plot2dModule:
					self.updateAxis(ax[2][0], title, [modVals[0], modVals[1]] if plotAnalitical else  vals[...,1])
			self.updateProjAxis(ax[0], title, vals, 0, dt)
			self.updateProjAxis(ax[1], title, vals, 1, dt)
			if plot2dModule:
				self.updateProjAxis(ax[2], title, modVals, 0, dt)

			
		
	def afterUpdateValues(self, newTime):
		timeTitle = "Time %4.4f" % newTime
		numFig = 0
		for fig in self.figures:
			fig.canvas.draw()
			if(hasattr(fig, "timeText")):
				fig.timeText.set_text(timeTitle)
			else:
				fig.timeText = fig.text(.1,.03,timeTitle)
			if saveImages:
				#make name compatible with ffmpeg
				#ffmpeg -r 1 -i img%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
				#the above does not work
				#ffmpeg -r 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
				#convert HANGS!!
				#convert -antialias -delay 1x2 *.png mymovie.mp4
				imgname = "%4.4f" % newTime
				if(len(imgname) == 6):
					imgname = "0"+imgname
				imgname = imgname.replace(".", "")	
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img%s.png"%imgname))
			numFig +=1
		#import time
		#time.sleep(5)

	def finish(self):
		pass



	def fftplot(self, ax, vals):
		from common import getDz1, getDz0, nint
		
		numPoints = (nint[0]+2) * (nint[1]+2)
		Y=fft(vals)/(numPoints)
		plotVals = np.absolute(Y)
		
		F0=fftfreq(numPoints, getDz0())
		F1=fftfreq(numPoints, getDz1())
		plotX , plotY = np.meshgrid(F0, F1)
		ax.grid(True)
		ax.plot_wireframe(plotX, plotY , plotVals)
#		print("plotVals")
#		np.set_printoptions(threshold='nan')
#		print(plotVals)
#		ax.imshow(plotVals)			


	def addFFTAxis(self, title, vals):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection="3d")
		ax.set_title(title)
		self.axes[title] = ax
		self.fftplot(ax, vals)
		self.figures.append(fig)
		fig.subplots_adjust(right=0.8)
		plt.draw()
		plt.show(block=False)
		
	def updateFFTAxis(self, title, vals):
		ax = self.axes[title]
		ax.cla()
		ax.set_title(title)
		self.fftplot(ax,vals)




