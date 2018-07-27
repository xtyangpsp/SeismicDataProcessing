USAGE: fmfocal.py fmdatafile outputfile [-S|s]
Plot first motion picks (can be picked in dbloc2 for Antelope database) and let the user pick nodal planes by visually fitting the data. The user needs to pick two planes for each event. The program will ask, when both planes are picked, whether to save the result or redo the picking of nodal planes. The results of picked nodal planes will be saved to [outputfile], which can be plotted using psmeca in GMT.

First motion data file format:
	evid lon lat depth magnitude event-to-station-azimuth dip/incidence fm

First motion conventions:
	up/compressional - cc;
	down/dilational - dd;
	null - nn; 
	up close to null - cn;
	down close to null - dn;

Note: In Anatelope database table 'predarr', dip angle is positive downward from horizontal.

Credits:
	Yinzhi Wang (original framework for plotting first motion data)
	Xiaotao Yang (improvement of workflow and GUI usability for looping through events)

Contacts:
		Xiaotao Yang (stcyang@gmail.com)

####### Library dependency #####
This program needs the following Python libraries (need to be installed):
	mpl_toolkits.basemap
	numpy
	matplotlib
	scipy

##### Test data:
After installing the required libraries, add it to searching path for use at any locations. In the current folder, there is a demo data set: firstmotion_demodata.d. This can be used to test the code. Run in terminal
	$fmfocal.py firstmotion_demodata.d testoutput.d
