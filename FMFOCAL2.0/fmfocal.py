#!/usr/bin/env python
############################################
# USAGE: fmfocal.py fmdatafile outputfile [-S|s]
# Plot first motion picks and let the user pick nodal planes by visually fitting the data.
# The user needs to pick two planes for each event. The program will ask, when both planes
# are picked, whether to save the result or redo the picking of nodal planes. The results
# of picked nodal planes will be saved to [outputfile], which can be plotted using psmeca
# in GMT.
#
# First motion data file format: evid lon lat depth magnitude event-to-station-azimuth dip/incidence fm
# First motion conventions:
# 		up/compressional - cc;
#		down/dilational - dd;
#		null - nn; 
#		up close to null - cn;
#		down close to null - dn;
# Note: In Anatelope database table 'predarr', dip angle is positive downward from horizontal.
# Credits:
#		Yinzhi Wang (original framework for plotting first motion data)
#		Xiaotao Yang (improvement of workflow and GUI usability for looping through events)
#
# Contacts:
# 		Xiaotao Yang (stcyang@gmail.com)
############################################
#
# ***** Code history ****
#Original version came from Yinzhi Wang.
#Modified by Xiaotao Yang started from May 7, 2013. Following is the list of major changes.
#
#Modified by Yinzhi Wang started from May 9, 2013. Following is the list of major changes.
#Data is read from file hardcoded in the program with the format as:
#EventID	Azimuth	Incidence	FirstMotion
#
#This program is now driven by the change of EventID
#
#Map projection changed to Azimuthal Equidistant Projection
#
#Turned plotting and mouse controlling into a function
#May 10, 2013      Xiaotao Yang     Print out the sequence number of event being processed.
#May 10, 2013      Yinzhi Wang     Changed color
#May 18-19, 2013	Xiaotao Yang	Read datafile from command line arguments. DoneChoice() to promot a window for user after each event.
#May 20-21, 2013	Xiaotao Yang	Read output file from command line arguments.
#					Write out the results into the output file for each event.
#					Changed first motion symbols/vonventions, see context for details.
#May 22, 2013		Xiaotao Yang	Changed input format to incude earthquake parameters. Defined class earthquakes to store the parameters.
#					Moved usage() and version() to the beginning to avoid loading libraries in case of usage error.
#Oct 18, 2013		Xiaotao Yang
#					1. here the incidence angle is equal to the dip angle. fixed the bug. 
#						dip angle in predarr of antelope database is positive downward from horizontal.
#					2. add P & T to output file.
#Oct 21, 2013		Xiaotao Yang
#					1. fixed a bug that skipped the last value for each event.
#Oct 31, 2013		Xiaotao Yang
#					1. added option to plot stations on basemap. -S|s
#July 27, 2018		Xiaotao Yang
# 					1. fixed bug error when calling m.drawparallels() by change dashes=[1,0] to dashes=[1,1];

import sys
global plotstation

def usage():
	print('** USAGE: fmfocal.py fmdatafile outputfile [-S|s]')
	print('** 		-S|s: plot station names')
	print('** Version: %s'%(version()))
def version():
	v='2.2'
	return v

if len(sys.argv)<3:
	usage()
	sys.exit(0)
elif len(sys.argv)==3:
	plotstation=False
elif len(sys.argv)==4:
	if sys.argv[3] == '-S' or sys.argv[3] =='-s':
		plotstation=True
else:
	print('**Wrong argument!')
	usage()
	sys.exit(0)

from mpl_toolkits.basemap import Basemap
import Tkinter
import tkMessageBox
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch
import math

class MouseMonitor(object):
	def __init__(self,quakes):
		self.quakes=quakes

	event = None
	lonlist = []
	latlist = []
	def mycall(self, event):
		self.event = event
		lonpt, latpt = m(event.xdata,event.ydata,inverse=True)
		lonpt=round(lonpt)
		latpt=round(latpt)
		self.lonlist.append(lonpt)
		self.latlist.append(latpt)
		#print 'LON = %s and LAT = %s' % (lonpt,-latpt)
		#ax = gca()  # get current axis
		#ax.hold(True) # overlay plots.
		# Plot a red circle where you clicked.
		#if(len(self.lonlist)<=4):
		x, y = m(lonpt,latpt)
		#ax.plot(x,y,'k.')
		m.scatter(x,y,11,marker='.',color='k')
		draw()  # to refresh the plot.
                #print len(self.lonlist)

		if(len(self.lonlist)%5==2):
			decl = calculate_declination(self.lonlist[-1],self.latlist[-1],self.lonlist[-2],self.latlist[-2])
                        #great(m, self.lonlist[-1], self.latlist[-1], decl, color='blue')
			great_circle_path(m, self.lonlist[-1],self.latlist[-1],self.lonlist[-2],self.latlist[-2], color='black')
			if decl<180:
				#great(m, self.lonlist[-1], self.latlist[-1], decl+180, color='blue')
				lonpt,latpt, backazimuth = shoot(self.lonlist[-1], self.latlist[-1], (decl+90+360)%360, math.pi*6378.13*0.5)
				x, y = m(lonpt,latpt)			
				print ('P1STRIKE = %s and P1PLUNGE = %s' % ((lonpt+360)%360,-latpt))
				self.lonlist.append(lonpt)
				self.latlist.append(latpt)
				m.scatter(x,y,13,marker='^',color='black')
			else:
				#great(m, self.lonlist[-1], self.latlist[-1], decl-180, color='blue')
				lonpt,latpt, backazimuth = shoot(self.lonlist[-1], self.latlist[-1], (decl-90+360)%360, math.pi*6378.13*0.5)
				x, y = m(lonpt,latpt)			
				print ('P1STRIKE = %s and P1PLUNGE = %s' % ((lonpt+360)%360,-latpt))
				self.lonlist.append(lonpt)
				self.latlist.append(latpt)
				m.scatter(x,y,13,marker='^',color='black')
			draw()
		elif(len(self.lonlist)%5==4):
			decl = calculate_declination(self.lonlist[-1],self.latlist[-1],self.lonlist[-2],self.latlist[-2])
			#great(m, self.lonlist[-1], self.latlist[-1], decl, color='blue')
			great_circle_path(m, self.lonlist[-1],self.latlist[-1],self.lonlist[-2],self.latlist[-2], color='black')
			if decl<180:
				#great(m, self.lonlist[-1], self.latlist[-1], decl+180, color='blue')
				lonpt,latpt, backazimuth = shoot(self.lonlist[-1], self.latlist[-1], (decl+90+360)%360, math.pi*6378.13*0.5)
				x, y = m(lonpt,latpt)			
				print ('P2STRIKE = %s and P2PLUNGE = %s' % ((lonpt+360)%360,-latpt))
				self.lonlist.append(lonpt)
				self.latlist.append(latpt)
				m.scatter(x,y,13,marker='s',color='black')
			else:
				#great(m, self.lonlist[-1], self.latlist[-1], decl-180, color='blue')
				lonpt,latpt, backazimuth = shoot(self.lonlist[-1], self.latlist[-1], (decl-90+360)%360, math.pi*6378.13*0.5)
				x, y = m(lonpt,latpt)			
				print ('P2STRIKE = %s and P2PLUNGE = %s' % ((lonpt+360)%360,-latpt))
				self.lonlist.append(lonpt)
				self.latlist.append(latpt)
				m.scatter(x,y,13,marker='s',color='black')
			lonpt1,latpt1, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], 90, math.pi*6378.13*0.5)
			lonpt2,latpt2, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], 180, math.pi*6378.13*0.5)
			lonpt3,latpt3, backazimuth = shoot(self.lonlist[-1], self.latlist[-1], 90, math.pi*6378.13*0.5)
			lonpt4,latpt4, backazimuth = shoot(self.lonlist[-1], self.latlist[-1], 180, math.pi*6378.13*0.5)
			backazimuth = calculate_arc(lonpt1, latpt1, self.lonlist[-1], self.latlist[-1])
			backazimuth1 = calculate_arc(lonpt3, latpt3, self.lonlist[-3], self.latlist[-3])
			plstrike1=(lonpt1+360)%360
			plstrike2=(lonpt3+360)%360
			pldip1=-latpt2
			pldip2=-latpt4
			plrake1=-backazimuth
			plrake2=-backazimuth1
			
			print ('PL1STRIKE = %s and PL1DIP = %s and PL1RAKE = %s' % (plstrike1,pldip1,plrake1))
			print ('PL2STRIKE = %s and PL2DIP = %s and PL2RAKE = %s' % (plstrike2,pldip2,plrake2))
			draw()

			decl = calculate_declination(self.lonlist[-3], self.latlist[-3], self.lonlist[-1], self.latlist[-1])
			lonpt,latpt, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], (decl+180)%360, math.pi*6378.13*0.25)
			x, y = m(lonpt,latpt)
			tstrike1=(lonpt+360)%360
			tplunge1=-latpt
			print ('TSTRIKE = %s and TPLUNGE = %s' % (tstrike1,tplunge1))
			#ax.plot(x,y,'yo')

			lonpt,latpt, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], decl, math.pi*6378.13*0.75)
			x, y = m(lonpt,latpt)
			tstrike2=(lonpt+360)%360
			tplunge2=-latpt
			print ('TSTRIKE = %s and TPLUNGE = %s' % (tstrike2,tplunge2))
			#ax.plot(x,y,'yo')

			lonpt,latpt, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], decl, math.pi*6378.13*0.25)
			x, y = m(lonpt,latpt)
			pstrike1=(lonpt+360)%360
			pplunge1=-latpt
			print ('PSTRIKE = %s and PPLUNGE = %s' % (pstrike1,pplunge1))
			#ax.plot(x,y,'ro')

			lonpt,latpt, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], (decl+180)%360, math.pi*6378.13*0.75)
			x, y = m(lonpt,latpt)
			pstrike2=(lonpt+360)%360
			pplunge2=-latpt
			print ('PSTRIKE = %s and PPLUNGE = %s' % (pstrike1,pplunge1))
			#ax.plot(x,y,'ro')
			
			if decl<180:
				lonpt,latpt, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], (decl+90+360)%360, math.pi*6378.13*0.5)
				x, y = m(lonpt,latpt)			
				#ax.plot(x,y,'go')
			else:
				lonpt,latpt, backazimuth = shoot(self.lonlist[-3], self.latlist[-3], (decl-90+360)%360, math.pi*6378.13*0.5)
				x, y = m(lonpt,latpt)			
				#ax.plot(x,y,'go')

			bstrike=(lonpt+360)%360
			bplunge=-latpt
			print ('BSTRIKE = %s and BPLUNGE = %s' % (bstrike,bplunge))
			#draw()
			
			draw()  # to refresh the plot.
		#else:
			# Close figure and go next or redo the current event. This is a GUI dialog window.
			#DoneChoice()
			result=tkMessageBox.askquestion("Satisfy?", "Save, close the figure and go next? Click NO to redo.",icon='warning')

			if result=='yes':
				#OUTPUTFILE HEADER: EVID	STRIKE	DIP	RAKE	TSTRIKE	TPLUNGE	PSTRIKE	PPLUNGE	BSTRIKE	BPLUNGE
#				fo.write('%d	%.4f	%.4f	%.4f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n' \
#					% (self.quakes.evid,self.quakes.lon,self.quakes.lat,self.quakes.depth,self.quakes.magnitude,\
#					plstrike1,pldip1,plrake1,tstrike1,tplunge1,pstrike1,pplunge1,bstrike,bplunge))
				fo.write('%d	%.4f	%.4f	%.4f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n' \
					% (self.quakes.evid,self.quakes.lon,self.quakes.lat,self.quakes.depth,self.quakes.magnitude,\
					plstrike1,pldip1,plrake1,plstrike2,pldip2,plrake2,pstrike1,pplunge1,tstrike1,tplunge1))
				plt.savefig("evid%d.ps"%(self.quakes.evid))
				close()
		#
			elif result=='no':
				print('Redo ...')
#end of class MouseMonitor.

class earthquakes:
	evid = []
	lon = []
	lat = []
	depth = []
	magnitude = []
# end of class earthquakes

			
def calculate_declination(lng_a, lat_a, lng_b, lat_b):
	
	degrees_to_radians = math.pi/180.0

	lat_a=lat_a*degrees_to_radians
	lng_a=lng_a*degrees_to_radians
	lat_b=lat_b*degrees_to_radians
	lng_b=lng_b*degrees_to_radians

	d=math.atan2(math.sin(lng_b-lng_a) * math.cos(lat_b),math.cos(lat_a)*math.sin(lat_b)-math.sin(lat_a)*math.cos(lat_b)*math.cos(lng_b-lng_a))/degrees_to_radians
	d=(d+360)%360
	return d

def calculate_arc(lng_a, lat_a, lng_b, lat_b):
	degrees_to_radians = math.pi/180.0

	phi1 = (90.0 - lat_a)*degrees_to_radians
	phi2 = (90.0 - lat_b)*degrees_to_radians

	theta1 = lng_a*degrees_to_radians
	theta2 = lng_b*degrees_to_radians
	
	arc= math.acos(math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + math.cos(phi1)*math.cos(phi2))
	
	#print 'Distance =	',arc*6371
	#print 'Colatitude =	',arc/degrees_to_radians
	
	return arc/degrees_to_radians

def shoot(lon, lat, azimuth, maxdist=None):
	#"""Shooter Function
	#Original javascript on http://williams.best.vwh.net/gccalc.htm
	#Translated to python by Thomas Lecocq
	#"""
	glat1 = lat * np.pi / 180.
	glon1 = lon * np.pi / 180.
	s = maxdist / 1.852
	faz = azimuth * np.pi / 180.
 
	EPS= 0.00000000005
	if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
		alert("Only N-S courses are meaningful, starting at a pole!")
 
	a=6378.13/1.852
	f=1/298.257223563
	r = 1 - f
	tu = r * np.tan(glat1)
	sf = np.sin(faz)
	cf = np.cos(faz)
	if (cf==0):
		b=0.
	else:
		b=2. * np.arctan2 (tu, cf)
 
	cu = 1. / np.sqrt(1 + tu * tu)
	su = tu * cu
	sa = cu * sf
	c2a = 1 - sa * sa
	x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
	x = (x - 2.) / x
	c = 1. - x
	c = (x * x / 4. + 1.) / c
	d = (0.375 * x * x - 1.) * x
	tu = s / (r * a * c)
	y = tu
	c = y + 1
	while (np.abs (y - c) > EPS):
 
		sy = np.sin(y)
		cy = np.cos(y)
		cz = np.cos(b + y)
		e = 2. * cz * cz - 1.
		c = y
		x = e * cy
		y = e + e - 1.
		y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
			  d / 4. - cz) * sy * d + tu
 
	b = cu * cy * cf - su * sy
	c = r * np.sqrt(sa * sa + b * b)
	d = su * cy + cu * sy * cf
	glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
	c = cu * cy - su * sy * cf
	x = np.arctan2(sy * sf, c)
	c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
	d = ((e * cy * c + cz) * sy * c + y) * sa
	glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi	
 
	baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
	glon2 *= 180./np.pi
	glat2 *= 180./np.pi
	baz *= 180./np.pi
 
	return (glon2, glat2, baz)

def great_circle_path(m, lng_a, lat_a, lng_b, lat_b, *args, **kwargs):
	degrees_to_radians = math.pi/180.0
	lat_a=lat_a*degrees_to_radians
	lng_a=lng_a*degrees_to_radians
	lat_b=lat_b*degrees_to_radians
	lng_b=lng_b*degrees_to_radians

	P=np.array([math.cos(lng_a)*math.cos(lat_a),math.sin(lng_a)*math.cos(lat_a),math.sin(lat_a)])
	Q=np.array([math.cos(lng_b)*math.cos(lat_b),math.sin(lng_b)*math.cos(lat_b),math.sin(lat_b)])
	w=Q-np.dot(P,Q)*P
	Qp=w/math.sqrt(np.dot(w,w))

	#print P
	#print Q
	#print w
	#print Qp

	latlist = []
	lnglist = []
	for t in range(362):
		g=math.cos(t*degrees_to_radians)*P+math.sin(t*degrees_to_radians)*Qp
		lng=math.atan2(g[1],g[0])/degrees_to_radians
		#lat=math.acos(g[2]/math.sqrt(np.dot(g,g)))/degrees_to_radians
		lat=math.atan2(math.sqrt(g[0]*g[0]+g[1]*g[1]),g[2])/degrees_to_radians
		lnglist.append(lng)
		latlist.append(90-lat)

		#g=cos(math.pi/-t)*P+sin(math.pi/-t)*Qp
		#lat=math.atan(g[1]/g[0])/degrees_to_radians
		#lng=math.acos(g[2]/math.sqrt(np.dot(g,g)))/degrees_to_radians
		#lnglist.append(lng)
		#latlist.append(90-lat)
	x=[]
	y=[]
	for i in range(len(lnglist)-1):
		xt, yt = m(lnglist[i],latlist[i])
		x.append(xt)
		y.append(yt)
		#m.scatter(xt,yt,size,marker='o',color='g')
		#print lnglist[i],latlist[i]
	ax = gca()  # get current axis
	ax.hold(True) # overlay plots.
	# Plot a red circle where you clicked.
	ax.plot(x,y, **kwargs)
	draw()
	return

# function to plot each event and call the mouse monitor object.
# first motion conventions: up/compressional - cc; down/dilational - dd; null - nn; up close to null - cn; down close to null - dn;
def plot_event(azimuth,incidence,motion,station, quakes):
        global m
        # create new figure, axes instances.
        fig=plt.figure(figsize=(6,6))
#        ax  = fig.add_subplot(1, 1, 1)
#	for child in ax.get_children():
#    		if isinstance(child, matplotlib.spines.Spine):
#        		child.set_color('#dddddd')
#        ax = plt.subplot(111)
        #ax1 = fig.add_axes([0.5,0.5,.5,.5],polar=True)
#        circle = Circle((0, 0), 4, facecolor='none', edgecolor=(0,0,0), linewidth=3, alpha=0.5)
#	ax.add_patch(circle)
#	ax.add_patch(circle)
	#plt.axis('off')
	#fig.add_subplot(111,projection='polar')
        plt.subplots_adjust(top=0.85)
        plt.title("Event ID: %d\n\n"%(quakes.evid))
	#plt.text(0.5, 1.08, "Event ID: %s"%(evid),horizontalalignment='center',fontsize=20)
        
        # setup map projection.
        #m = Basemap(projection='robin',lon_0=0)
        #m = Basemap(projection='spstere',boundinglat=0,lon_0=180,resolution='l')
        #m = Basemap(width=20000000,height=20000000,projection='aeqd',lat_0=-90,lon_0=0)
        m = Basemap(projection='spaeqd',boundinglat=0,lon_0=180)
        # draw parallels
        #m.drawparallels(np.arange(-90,10,10),color='0.75')
        m.drawparallels([0,0],color='0',dashes=[1,1])

        # draw meridians
        #m.drawmeridians(np.arange(0,360,10),color='0.75')
#	im=plt.show()
#        im.set_clip_path(circle)
        
        size=50
        textfontsize=8
        for i in range(len(azimuth)):
                x, y = m(azimuth[i],-(incidence[i]))
                #print('incidence:%f' %(-(incidence[i])))
                #print('%f' %y)
                #plt.text(x,y,station[i],fontsize=9, ha='center',va='center',color='r')
                if (motion[i]=='cc'):
                        m.scatter(x,y,size,marker='o',color='k')
                        if(plotstation):
                        	plt.text(x,y+500000,station[i],fontsize=textfontsize, ha='center',va='center',color='k')
                elif (motion[i]=='dd'):
                        m.scatter(x,y,size,marker='o',facecolors='none', edgecolors='k')
                        if(plotstation):
                        	plt.text(x,y+500000,station[i],fontsize=textfontsize, ha='center',va='center',color='r')
                elif (motion[i]=='nn'):
                        m.scatter(x,y,size,marker='x',color='k')
                        if(plotstation):
                        	plt.text(x,y+500000,station[i],fontsize=textfontsize, ha='center',va='center',color='b')
                elif (motion[i]=='cn'):
                        m.scatter(x,y,size*0.7,marker='o',color='k',)
                        m.scatter(x,y,size,marker='x',color='k')
                elif (motion[i]=='dn'):
                        m.scatter(x,y,size,marker='o',facecolors='none', edgecolors='k')
                        m.scatter(x,y,size,marker='x',color='k')
        mouse = MouseMonitor(quakes)
        connect('button_press_event', mouse.mycall)
        plt.show()

def DoneChoice():
#Close figure and go next or redo the current event.
# Currently, this funciton is embeded in MouseMonitor() for easy loop through events.
	result=tkMessageBox.askquestion("Go Next?", "Do you want to close the figure and go next? Click NO to redo.",icon='warning')

	if result=='yes':
		close()
		#
	elif result=='no':
		print('Redo ...')

#-------------------------------------------Main function------------------------------------------#
# open first motion data file: evid lon lat depth magnitude event-to-station-azimuth dip/incidence fm
#first motion conventions: up/compressional - cc; down/dilational - dd; null - nn; up close to null - cn; down close to null - dn;
f=open(sys.argv[1])
fout=sys.argv[2]

global fo

#create output file object fo.
fo=open(fout,'w')
#fo.write('EVID	LON	LAT	DEPTH	MAG	STRIKE	DIP	RAKE	TSTRIKE	TPLUNGE	PSTRIKE	PPLUNGE	BSTRIKE	BPLUNGE\n')
fo.write('EVID	LON	LAT	DEPTH	MAG	STRIKE1	DIP1	RAKE1  STRIKE2	DIP2 RAKE2	PSTRIKE	PPLUNGE	TSTRIKE TPLUNGE\n')
data=f.readlines()
f.close()
data=[x.split() for x in data]
for j in range(len(data)):
        length=len(data)
        for i in range(len(data)):
                if data[i][0][0]=='#':
                        temp=data[0:i]
                        if temp==[]:
                                data=data[i+1:len(data)]
                        elif i+1!=len(data):
                                data=temp+data[i+1:len(data)]
                        else:
                                data=temp;
                        break
        length2=len(data)
        if length2==length:
                break
data=[[int(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4]),float(x[5]),float(x[6]),str(x[7]), str(x[8])] for x in data]
#print('length of data: %d'%(len(data)))
azimuth=[]
incidence=[]
motion=[]
station=[]
quakes=earthquakes
start=0
enum=0
for j in range(len(data)-1):
        if data[j][0]!=data[j+1][0]:
        	#print('line= %d' %(j))
        	#print('%d,  %d' %(data[j][0],data[j+1][0]))
                for i in range(start,j+1):
                	#print('i= %d' %i)
                        azimuth.append(data[i][5])
                        incidence.append(data[i][6])
                        motion.append(data[i][7])
                        station.append(data[i][8])
                enum=enum+1
                quakes.evid=data[j][0]
                quakes.lon=data[j][1]
                quakes.lat=data[j][2]
                quakes.depth=data[j][3]
                quakes.magnitude=data[j][4]
                print('DataforEventID:	%d			--> %d' %(quakes.evid,enum))
                
                plot_event(azimuth,incidence,motion,station,quakes)
                start=j+1
                azimuth=[]
                incidence=[]
                motion=[]
        elif j+2==len(data):
                for i in range(start,len(data)):
                	#print('i= %d' %i)
                        azimuth.append(data[i][5])
                        #print('%s' %(data[i][7]))
                        incidence.append(data[i][6])
                        motion.append(data[i][7])
                        station.append(data[i][8])
                enum=enum+1
                #print('%d' %len(azimuth))
                quakes.evid=data[j][0]
                quakes.lon=data[j][1]
                quakes.lat=data[j][2]
                quakes.depth=data[j][3]
                quakes.magnitude=data[j][4]
                print('DataforEventID:	%d			--> %d' %(quakes.evid,enum))
                plot_event(azimuth,incidence,motion,station,quakes)

fo.close()
