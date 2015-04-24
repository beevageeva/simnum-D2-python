import matplotlib
useWindow = True
if useWindow:
  matplotlib.use('TkAgg')
else:
  matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from time import sleep


def combineImage(imgname,trajname):
	from matplotlib.image import imread
	imgVals = imread(imgname)
	trajVals = imread(trajname)	
	masked_data = np.ma.masked_where(trajVals==1 , trajVals)
	plt.imshow(imgVals)
	plt.imshow(masked_data, cmap=cm.gray)
	plt.draw()
	plt.show(block=False)
	sleep(10)		

def combineAll():
	import glob
	for imgname in glob.glob('img*.png'):
		combineImage(imgname,"Trajectory.png")

#test
combineImage("img007213.png", "Trajectory.png")
#combineAll()
