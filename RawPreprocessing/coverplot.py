#!/usr/bin/env python
'''
Graphical tool to plot seismic waveform, and dataless coverage
'''
VERSION = '2013.295'

from Tkinter import *
import itertools
#CONSTANTS
WINDOW_WIDTH = 1000
WINDOW_WIDTH = 2000
WINDOW_HEIGHT = 800
WINDOW_HEIGHT = 1200
LEFTFRAME_WIDTH = 200
BOTTOMFRAME_HEIGHT = 50

#from flowlib import synclib

#debug
import sys
import time
TEST_SYNC = "/Users/lloyd/work/test_data/X2.2009.2011.dmc.sync.mshd.cshd"
TEST_DL = "/Users/lloyd/work/test_data/X2.2010.264.X2.10.SAHKE_db.20102570910.dataless"
#TEST_SYNC = "/Volumes/flow/AUTO/EXPS/Z9.2010/DATA/SYNC/Z9.2010.2012,340.dmc.sync"
#TEST_DL = "/Volumes/flow/AUTO/EXPS/Z9.2010/DATA/DATALESS/Z9.12.SESAME_S71_db.20123260846.dataless"

class CoverPlot(Frame):
    #constants
    def __init__(self, master, ms_syncs, dl_syncs) :
        #superclass init
        Frame.__init__(self, master)
        self.master = master
        #self.ms_sum_list = ms_sum_list
        #self.dl_sum_list = dl_sum_list
        
        #Menubar do we need?
        #self.menubar = Menu(self.master)
        #self.filemenu = Menu(self.menubar, tearoff=0)
        #self.menubar.add_cascade(label='File', menu=self.filemenu)
        #self.master.config(menu=self.menubar)
        
        self.pack(fill=BOTH, expand=YES)
        self.leftFrame = Frame(self, width=LEFTFRAME_WIDTH, borderwidth=2, relief=RIDGE)
        #test button
        self.ButtonPlot = Button(self.leftFrame, text="Plot", command=self.make_plot)
        self.ButtonPlot.pack()
        self.ButtonSave = Button(self.leftFrame, text='Save', command = self.save)
        self.ButtonSave.pack()
        self.plotFrame = Frame(self, bg='blue', borderwidth=2,  relief=RIDGE)
        self.bottomFrame = Frame(self, height=BOTTOMFRAME_HEIGHT, borderwidth=2,  relief=RIDGE)
        self.bottomFrame.pack(side=BOTTOM, expand=NO, fill=X)
        self.leftFrame.pack( side=LEFT, expand=NO, fill=Y)
        self.plotFrame.pack(side=LEFT,fill=BOTH, expand=YES)
        
        #def mouse_wheel_down(event) :
            #''' Button 4
            #scroll up like new mac style dragging the canvas down'''
            #self.canvas.yview_scroll(-1,'units')
        #def mouse_wheel_up(event) :
            #''' Button 5
            #scroll down like new mac style dragging the canvas up'''
            #self.canvas.yview_scroll(1,'units')
        
       #plot canvas 
        self.canvas = self.PlotCanvas(self.plotFrame, zoom_func=self.make_plot)
        #self.canvas = Canvas(self.plotFrame, bg='black', yscrollcommand=self.scrollbar.set)
        self.canvas.pack(side=LEFT, expand=YES, fill=BOTH)
        #self.canvas.bind(sequence='<Button-4>', func=lambda e: self.canvas.yview_scroll(-1, 'units'))
        #self.canvas.bind(sequence='<Button-4>', func=mouse_wheel_down)
        #self.canvas.bind(sequence='<Button-5>', func=mouse_wheel_up)
       
        #May need to put this somewhere
        #self.canvas.config(scrollregion=(0,0,X,Y)
        #??? debug
        self.ms_syncs = ms_syncs
        self.dl_syncs = dl_syncs
        self.load_data()
        
    # PlotCavas class encompasses all and only the plottying things
    class PlotCanvas(Canvas):
        #subclass of canvas to add other info to
        #margins
        MLEFT = 10
        MRIGHT = 80
        MTOP = 10
        MBOTTOM = 10
        ROW_HEIGHT = 30
        WIDE_ROW = 15
        THIN_ROW = 5
        MIN_LINE_WIDTH = 1
        # X value plottable area left and right
        X0 = MLEFT +200 
        def __init__ (self, master, zoom_func):
            self.masterFrame = master
            self.scrollbar = Scrollbar(self.masterFrame)
            Canvas.__init__(self, self.masterFrame, bg='black', yscrollcommand=self.scrollbar.set)
            self.scrollbar.pack(side=RIGHT, fill=Y)
            self.scrollbar.config(command=self.yview)
            #scrolling
            self.bind(sequence='<Button-4>', func=lambda e: self.yview_scroll(-1, 'units'))
            self.bind(sequence='<Button-5>', func=lambda e: self.yview_scroll(1, 'units'))
            #initialize
            #Scale on top
            self.timecanvas = Canvas(self.masterFrame, bg='black', height=self.ROW_HEIGHT *2)
            self.timecanvas.pack(side=TOP, expand=False, fill=X)
            #dict of sncl : index (order of sncl in plot
            self.master_index = None
            self.earliest = None
            self.latest = None
            #zooming shift click
            self.bind(sequence='<Shift-Button-1>',  func=self.zoom)
            self.Zclick1 = None
            self.zoom_func = zoom_func
            #info clicking ctrl-click
            self.bind(sequence='<Control-Button-1>', func=self.display_time)
        def zoom(self, event) :
            if not self.Zclick1:
                #first click save, and draw a yellow line
                self.Zclick1 = self.canvasx(event.x)
                self.create_line(self.Zclick1, 0, self.Zclick1, self.canvas_height, fill='yellow', width=2)
            else:
                #second click ->zoom
                click1 = self.x2epoch(self.Zclick1)
                #start = ( self.Zclick1 - self.X0 ) / self.pixelsPERsecond + self.earliest
                click2 = self.x2epoch(self.canvasx(event.x))
                self.Zclick1 = None
                start = min(click1, click2)
                end = max(click1, click2)
                self.zoom_func(start, end) 
        def display_time(self, event) :
                # for ctrl-click shows vert bar and time
                x = self.canvasx(event.x)
                y = self.canvasy(event.y)
                self.create_line(x, 0, x, self.canvas_height, fill='white', width=2)
                display_time = Sync.epoch2string(self.x2epoch(x))
                self.create_text(x, y, text=display_time, fill='white', anchor='w')
                
                
        def x2epoch(self, x) :
            #returns the epoch time of an x value e.g. convert clicks to time
            return ( x -self.X0 ) / self.pixelsPERsecond + self.earliest
            
            
        def setup_canvas(self, start=None, end=None):
            # run this once before each plot  after a window resize etc.
            self.delete(ALL)
            self.canvas_width = self.winfo_width()
            #self.canvas_width = self.masterFrame.winfo_width()
            self.XN = self.canvas_width - self.MRIGHT
            self.Xlen = self.XN - self.X0
            #build list of all sncls from all sums
            #masterlist = itertools.chain(*(s for s in list_sumdicts))
            ##uniquify and sort list
            #masterlist = sorted({}.fromkeys(masterlist).keys())
            ## index of row for a given sncl to be plotted
            #self.sncl_index = dict(zip(masterlist, itertools.count()))
            
            #setup  scale
            if start is not None:
                self.earliest = start
            if end is not None:
                self.latest = end
            #else:
                ## ?? temp just the first sumdict
                #self.earliest = min( ( x[0][0] for x in list_sumdicts[0].values() ) )
                #self.latest = max( ( x[-1][-1] for x in list_sumdicts[0].values() ) )
            self.pixelsPERsecond = float(self.Xlen) / (self.latest - self.earliest)
            
            #setup time scale
            self.setup_time_axes()
            
        def fmt_jul_time(self, epoch):
            #returns a julian string with enough detail for plot zoom level
            return time.strftime(self.timeformat_jul, time.gmtime(epoch) )
        def fmt_greg_time(self, epoch):
            return time.strftime(self.timeformat_jul, time.gmtime(epoch) )
        def fmt_jul_date(self,epoch):
            return time.strftime("%Y : %j", time.gmtime(epoch) )
        def fmt_greg_date(self,epoch):
            return time.strftime("%m / %d", time.gmtime(epoch) )
        def fmt_time(self, epoch):
            return time.strftime("%H : %M : %S", time.gmtime(epoch) )
        
        def setup_time_axes(self):
            #shorthand
            tc = self.timecanvas
            #some constants
            MIN = 60
            HOUR = MIN*60
            DAY = HOUR*24
            canvas_width = tc.winfo_width()
            canvas_height = tc.winfo_height()
            barY = canvas_height - 5
            #date time formater depends on zoom scale
            #default (very zoomed out only show day
            self.timeformat_jul = "%Y:%j"
            #show hour
            if self.pixelsPERsecond > 1.0 / (HOUR):
                self.timeformat_jul = "%Y:%j:%H"
            #show min
            if self.pixelsPERsecond > 1.0 / ( MIN):
                self.timeformat_jul = "%Y:%j:%H:%M"
            #show seconds
            if self.pixelsPERsecond > 1.0 :
                self.timeformat_jul = "%Y:%j:%H:%M:%S"
            
            #clear the scale
            tc.delete(ALL)
            tc.create_text(self.MLEFT, self.MTOP, text='%0.2f days' % ((float(self.latest) - self.earliest) / DAY ), fill='white', anchor='w' )
            #Display pixels per second debug
            #tc.create_text(self.MLEFT, self.MTOP*3, text="%f pix/sec" % self.pixelsPERsecond, fill='white', anchor='w' )
            #draw horizontal line
            tc.create_line(self.X0, barY, self.XN, barY, fill='white', width=1)
            
            #test? show start and end time
            #start
            #Ylabel = self.MTOP + 15
            #startlabel = tc.create_text(self.X0, Ylabel, text=self.fmt_jul_time(self.earliest), fill='white', anchor='s')
            #end
            #tc.create_text(self.XN, Ylabel, text=self.fmt_jul_time(self.latest), fill='yellow', anchor='s')
            
            #real labels
            num_labels = self.Xlen / 150
            label_width = self.Xlen / num_labels
            for x in range(self.X0, self.XN, label_width):
                tc.create_text( x, self.MTOP +10 , text=self.fmt_jul_date(self.x2epoch(x)), fill='white', anchor='s')
                tc.create_text( x, self.MTOP + 22, text=self.fmt_greg_date(self.x2epoch(x)), fill='white', anchor='s')
                if self.pixelsPERsecond > 1.0 / HOUR:
                    tc.create_text( x, self.MTOP + 34, text=self.fmt_time(self.x2epoch(x)), fill='white', anchor='s')
                tc.create_line(x,barY, x, barY - 5, fill='white', width=1 )
            
            

            
        def update_after_plot(self):
            self.config(scrollregion=( (0,0) + self.bbox(ALL)[2:4] ) )
            self.canvas_height = self.bbox(ALL)[3]
        
        def plot_labels(self, color = 'cyan' ):
            '''Plots channel names on left from self.master_index)
            call setup_canvas first once before plotting and after changes i.e. window resize'''
            for sncl in self.sncl_index:
                Y = self.sncl_index[sncl] * self.ROW_HEIGHT + self.MTOP
                self.create_text(self.MLEFT, Y, text=':'.join(sncl), fill=color, anchor='w')
                #self.create_line(self.X0, Y, self.X0, Y, fill=color, width=15)
                #self.create_line(self.X0+2, Y, self.X0+3, Y, fill=color, width=15)
                #self.create_line(self.X0+5, Y, self.X0+7, Y, fill=color, width=15)
        
        def plot_single_sum_dict(self, sum_dict, width, color):
            '''Plots a single sum_dict on self'''
            for sncl in sum_dict:
                Y = self.sncl_index[sncl] * self.ROW_HEIGHT + self.MTOP
                for start, end in sum_dict[sncl]:
                    #is off canvas?
                    if end <= self.earliest or start >= self.latest:
                        continue
                    if start < self.earliest:
                        start = self.earliest
                    if end > self.latest:
                        end = self.latest
                    #X1 = max (1, (start - self.earliest) * self.pixelsPERsecond + self.X0 )
                    #X2 = max (1, (end - self.earliest) * self.pixelsPERsecond + self.X0 )
                    # was not plotting min size lines
                    X1 = (start - self.earliest) * self.pixelsPERsecond + self.X0
                    X2 = max ( X1+ self.MIN_LINE_WIDTH, (end - self.earliest) * self.pixelsPERsecond + self.X0 )
                    self.create_line(X1, Y, X2, Y, fill=color, width=width)
                    
    # End PlotCanvas
        
    def load_data(self):
        # loads from 2 synclists ms and dl
        #this is a list of tuples ( coverage_dict, gap_dict ) of epoch dicts
        earlylate = [sys.maxint, 0]
        self.ms_dicts = list() 
        for sync in self.ms_syncs :
            gd = dict()
            sd = sync.mk_sum_dict(1, gap_dict=gd, earlylate=earlylate)
            self.ms_dicts.append((sd,gd))
            
        self.dl_dicts = list()
        for sync in self.dl_syncs :
            gd = dict()
            sd = sync.mk_sum_dict(1, gap_dict=gd, earlylate = earlylate)
            self.dl_dicts.append((sd,gd))
        
        self.earliest, self.latest = earlylate
            
        masterlist = list()
        for l in self.ms_dicts, self.dl_dicts:
            for t in l:
                masterlist.extend( t[0].keys() )
        #list_sumdicts = [ s[0].keys() for s in( d for d in ( l for l in (self.dl_dicts, self.ms_dicts) ) ) if len(s) > 0 ]
        ##masterlist = itertools.chain(*)
        #masterlist = itertools.chain( * ( s[0] for s in ( d for d in ( l for l in (self.dl_dicts, self.ms_dicts) ) ) if len( s) > 0 ) )
        ##masterlist = itertools.chain(*(s for s in list_sumdicts))
        #uniquify and sort list
        masterlist = sorted({}.fromkeys(masterlist).keys())
        # index of row for a given sncl to be plotted
        self.sncl_index = dict(zip(masterlist, itertools.count()))
        self.canvas.sncl_index = self.sncl_index
        
    def make_plot(self, start=None, end=None):
        #button was pressed
        thin = 13
        thick = 17
        start = start or self.earliest
        end = end or self.latest
        self.canvas.setup_canvas(start=start, end=end)
        
        self.canvas.plot_labels()
        for sd in ( s[0] for s in self.ms_dicts):
            #wide data in red will only show when no dl plotted over
            self.canvas.plot_single_sum_dict(sd, thick, 'red')
        for dl_sd in (d[0] for d in self.dl_dicts):
            #wide DL plotted
            self.canvas.plot_single_sum_dict(dl_sd, thick, 'deep sky blue')
        for dl_gd in ( d[1] for d in self.dl_dicts):
            # black out DL gap to cover ms gap == no dl and no MS = black
            self.canvas.plot_single_sum_dict(dl_gd, thin, 'black')
        for sd in (s[0] for s in self.ms_dicts):
            # wf coverage
            self.canvas.plot_single_sum_dict(sd, thin, 'dark blue')
        for gd in (s[1] for s in self.ms_dicts):
            #waveform gap
            self.canvas.plot_single_sum_dict(gd, thin, 'deep pink')
        
        self.canvas.update_after_plot()
        
    def save( self, filename='plot.eps'):
        print "Saving plot to %s ... " % filename, 
        self.canvas.plot_labels(color='black')
        self.canvas.postscript(file=filename, width=self.canvas.canvas_width, height=self.canvas.canvas_height, x=0, y=0 )
        print "Complete."
        self.canvas.plot_labels()
        
    def save_png(self, filename='plot.png'):
        print "Saving plot to %s ... " % filename, 
            
       
    def plot_sum_dict(self) :
        '''Plots sum_dict from Sync on self.canvas'''
        #setup for CANVAS
        
        
        #?????????? this doesn't shrink after the first plot!!!! use the frame and subtract!!!
        self.canvas.delete(ALL)
        canvas_width = self.plotFrame.winfo_width()
        #canvas_width = self.canvas.winfo_width()
        #canvas_height = self.canvas.winfo_height()
        #margins
        MLEFT = 10
        MRIGHT = 10
        MTOP = 10
        MBOTTOM = 10
        ROW_HEIGHT = 30
        WIDE_ROW = 15
        THIN_ROW = 5
        # X value plottable area left and right
        X0 = MLEFT + 150
        XN = canvas_width - MRIGHT
        
        #----------test stuff------------------#
        #make a sync sum_dict
        from flowlib import synclib
        s = Sync()
        s.read_file(TEST_SYNC)
        gd = dict()
        sd = s.mk_sum_dict(1,gap_dict=gd)
        #lets dry a dl
        sdl = Sync()
        sdl.read_dl_file(TEST_DL)
        dl_gd = dict()
        dl_sd = sdl.mk_sum_dict(1, gap_dict=dl_gd)
        #build a master set of all keys
        master_keys = list()
        for sum in itertools.chain(sd, dl_sd):
            #master_keys.extend(sum)
            #these aren't lists of sum_dicts just sum_dicts
            master_keys.append(sum)
        master_keys = sorted({}.fromkeys(master_keys).keys())
        master_index = dict(zip(master_keys, range(len(master_keys))))
        #setup
        earliest = min( ( x[0][0] for x in sd.values() ) )
        latest = max( ( x[-1][-1] for x in sd.values() ) )
        pixelsPERsecond =float(XN - X0) / (latest - earliest)
        
        #plot dl
        for sncl in dl_sd.keys():
            Y = master_index[sncl] * ROW_HEIGHT + MTOP
            self.canvas.create_text(MLEFT, Y, text = ':'.join(sncl), fill='cyan', anchor='w')
            for start, end in  dl_sd[sncl]:
                X1 = (start - earliest) * pixelsPERsecond + X0
                X2 = max(X1+1, (end - earliest) * pixelsPERsecond + X0)
                #X2 = (end - earliest) * pixelsPERsecond + X0
                #print X1, X2
                #draw coverage line
                self.canvas.create_line(X1, Y, X2, Y, fill='blue', width=15 )
        #plot dl gaps
        for sncl in dl_gd.keys():
            Y = master_index[sncl] * ROW_HEIGHT + MTOP
            self.canvas.create_text(MLEFT, Y, text = ':'.join(sncl), fill='cyan', anchor='w')
            for start, end in  dl_gd[sncl]:
                X1 = (start - earliest) * pixelsPERsecond + X0
                X2 = max(X1+1, (end - earliest) * pixelsPERsecond + X0)
                #X2 = (end - earliest) * pixelsPERsecond + X0
                #print X1, X2
                #draw coverage line
                self.canvas.create_line(X1, Y, X2, Y, fill='orange', width=18 )
        #plot wf
        for sncl in sd.keys():
            Y = master_index[sncl] * ROW_HEIGHT + MTOP
            self.canvas.create_text(MLEFT, Y, text = ':'.join(sncl), fill='cyan', anchor='w')
            for start, end in  sd[sncl]:
                X1 = (start - earliest) * pixelsPERsecond + X0
                X2 = max(X1+1, (end - earliest) * pixelsPERsecond + X0)
                #X2 = (end - earliest) * pixelsPERsecond + X0
                #print X1, X2
                #draw coverage line
                self.canvas.create_line(X1, Y, X2, Y, fill='darkblue', width=5 )
        #plot gap
        for sncl in gd.keys():
            Y = master_index[sncl] * ROW_HEIGHT + MTOP
            self.canvas.create_text(MLEFT, Y, text = ':'.join(sncl), fill='cyan', anchor='w')
            for start, end in  gd[sncl]:
                X1 = (start - earliest) * pixelsPERsecond + X0
                X2 = max(X1+3, (end - earliest) * pixelsPERsecond + X0)
                #X2 = (end - earliest) * pixelsPERsecond + X0
                #print X1, X2
                #draw coverage line
                self.canvas.create_line(X1, Y, X2, Y, fill='red', width=7 )
            
        self.canvas.config(scrollregion=(self.canvas.bbox(ALL)))
        pass
        
                
def get_command_args():
    '''get ms_sum_list and dl_sum_list from command line'''
    import optparse
    parser = optparse.OptionParser(usage="Usage: %prog [options]", version=VERSION)
    parser.add_option('-m', '--ms', action='append',
                      help= 'MiniSEED sync file' )
    parser.add_option('-d', '--dl', action='append', 
                      help = 'DatalessSEED file. If given will be compared with MSEED sync given, and a summary of the MSEED sync will not be produced.' )
    (options, args) = parser.parse_args()
    print options
    print args
    dl_syncs = list()
    ms_syncs = list()
    for dl in options.dl or ():
        '''Dataless sync given: compare'''
        '''tolerance command line arg is not used for dataless, always use zero and print day gran'''
        print "Reading DL %s" % dl
        dl_sync = Sync()
        dl_sync.read_dl_file(dl)
        dl_syncs.append(dl_sync)

    for ms in options.ms or ():
        '''no dataless sync given: produce summary'''
        insync = Sync()
        try:
            ms_crshd = "%s.crshd" % ms
            print "Attempting to read %s ..." % ms_crshd , 
            insync.read_file(ms_crshd)
            print "Success"
        except:
            print
            print  "%s doesn't exist. Reading original sync %s" % (ms_crshd, ms)
            insync.read_file(ms)
            print "Crushing nonwf channels out of sync..."
            insync.crush_nonwf()
            print "Writing crushed sync %s" % ms_crshd
            insync.write_file( ms_crshd)
        ms_syncs.append(insync)

    return ms_syncs, dl_syncs
    
    
### END of coverplot.py

### BEGIN synclib.py

'''
synclib.py :: All things related to dmc sync files
references http://www.iris.washington.edu/SeismiQuery/goat/syncformat.html
--Lloyd Carothers
'''
VERSION = "2013.022"
import re
import os.path
import sys
from math import ceil
import itertools
from time import localtime, mktime, strptime, gmtime, strftime
from calendar import timegm
import subprocess

#Log crushing tolerance 1/2 a day
LOGTOLERANCE = 60*60*24/2
#define sync fields
NET = 0     #Network
STA = 1     #Station
LOC = 2     #Location
CHAN = 3    #Channel
ST = 4      #Start time
ET = 5      #End time
MCD = 6     #Max clock drift
SR = 7      #Sample rate
NS = 8      #Number of samples
CF = 9      #Channel flag
SV = 10     #Station volume
TN = 11     #DCC tape number
VN = 12     #DMC volume number
COM = 13    #Comment
DLMDMC = 14 #Date line modified by DMC
DLMDCC = 15 #Date line modified by DCC
#define time fields
YEAR = 0
DAY = 1
HOUR = 2
MIN = 3
SEC = 4

#indexes for keys
kNET = 0
kSTA = 1
kLOC = 2
kCHAN = 3
kSR = 4

#indexes for keys with date
kdNET = 0
kdSTA = 1
kdLOC = 2
kdCHAN = 3
kdST = 4
kdSR = 5

class Sync (object):
    '''
    syncfile object
    can be read, altered and writen out
    subclassed list:
        list of tuplesc

    dl_times
    functions that end in _dlt are for using dl_times output to create a sync output
    dl_times has not sample rate so comparisons and keys don't use it
    '''
    STP = '%Y,%j,%H:%M:%S'

    def __init__(self, string=None):
        #list of self is a list of tuples of sync fields
        self.slist = None
        #can be called with a string that will populate self
        if string is not None:
            self.read_string(string)
        # epochdict can be made once and reused
        self.epochdict = None
        self.epochdict_dlt = None

    def __eq__(self, other):
        #are the types the same
        if type(self) is type(other):
            #are the lengths of the list the same
            if len(self. slist) != len(other.slist):
                return False
            #is each line the same
            for index, line in enumerate(self.slist):
                if line[:6] != other.slist[index][:6]:
                    return False
            return True

        else:
            return False

    def show_head(self):
        for line in self.slist[:30]:
            print line
    def populate(self, iterable):
        #iteralble should give a string line on .next()
        #print "\tReading"
        self.slist = [ t for t in (line.strip().split('|') for line in iterable) if len(t)>=16 ]
    def populate_dlt(self, iterable, net):
        '''Populate from dl_times output, NOT a DatalessSEED volume'''
        self.slist = list()
        sta = None
        #line is either STA: or {wspace}CHAN:LOC {wspace} start_time {wspace} end_time
        for line in iterable:
            if len(line) < 1 :continue
            assert isinstance(line, (str, unicode))
            if not line[0].isspace():
                #station line
                #remove colon
                sta = line.strip()[:-1]
                continue
            else:
                #channel line
                if sta is None:
                    raise RuntimeError, "Channel line has no station %s" % line
                fields = line.strip().split()
                #gives (CHAN:LOC, start, end)
                #channel and loc
                #gives (chan, loc) loc may be ''
                chan, loc = fields[0].split(':')
                #start and end times, remove fraction of second
                st = fields[1].split('.')[0]
                et = fields[2].split('.')[0]
                self.slist.append( [net, sta, loc, chan, st, et] )

    def sort_mash(self):
        #print "\tSorting"
        self.slist.sort()
        #print "\tRemoving dups"
        tot = len(self.slist)
        self.remove_duplicates()
        #print "\tRemoved ", tot - len(self.slist)
        #print "\tMashing non_wf"
        tot = len(self.slist)
        self.mash_nonwf()
        #print "\tRemoved", tot - len(self.slist)
        #print "\tLines Remaining", len(self.slist)
    def sort_crush(self):
        #print "\tSorting"
        self.slist.sort()
        #print "\tRemoving dups"
        tot = len(self.slist)
        self.remove_duplicates()
        #print "\tRemoved ", tot - len(self.slist)
        #print "\tMashing non_wf"
        tot = len(self.slist)
        self.crush_nonwf()
        #print "\tRemoved", tot - len(self.slist)
        #print "\tLines Remaining", len(self.slist)
    def sort_mash_keep_dups(self):
        self.slist.sort()
        self.mash_nonwf()

    #READERS
    def read_file(self, filename):
        #print "Using", filename
        infh = open(filename, 'r')
        self.populate(infh)
        infh.close()
    def read_dl_file(self, filename):
        '''Read a datalessSEED volume NOT a dl_times output
        This means we can get the sample rate and network code, 2 things which dl_times drops'''
        #These are the SEED packet that print in the beginning of rdseeds output grepping boost performance significantly
        PACKETS = ( "B050F03", "B052F03", "B052F04", "B052F22", "B052F23", "B052F18", "B050F16" )
        rdseed_args = ["rdseed", '-s', '-f', filename]
        grep_args = ["grep", "-E", '|'.join(PACKETS)]
        #print os.environ['PATH']
        try:
            process1 = subprocess.Popen(rdseed_args, stdout=subprocess.PIPE)
            process2 = subprocess.Popen(grep_args, stdin=process1.stdout, stdout=subprocess.PIPE)
        except OSError, e:
            sys.stderr.write("Error: Could not run rdseed\n%s\n" % e)
            sys.exit(-1)
        out, err = process2.communicate()
        if process2.returncode != 0:
            sys.stderr.write("Could not run rdseed: %s\n" % args)
            sys.exit(-1)
        self.slist = list()
        #sentinel for channels that have no end times
        null_end_time = False
        for line in out.splitlines():
            if re.match('B050F03', line):
                #Station Code:
                sta = line.split()[3]
                continue
            if re.match('B050F16', line):
                #Network Code
                net = line.split()[3]
                continue
            if re.match('B052F04', line):
                #Channel:
                chan = line.split()[2]
                continue
            if re.match('B052F03', line):
                #Location:
                words = line.split()
                if len(words) >2:
                    loc = line.split()[2]
                else:
                    loc = ''
                continue
            if re.match('B052F18', line):
                #Sample rate:
                sr = line.split()[3]
                continue
            if re.match('B052F22', line):
                #Start date:
                start = line.split()[3]
                continue
            if re.match('B052F23', line):
                #End date:
                end = line.split()[3]
                if end == "(null)":
                    sys.stderr.write("No end time for %(net)s.%(sta)s.%(loc)s.%(chan)s.%(start)s\n" % locals() )
                    null_end_time = True
                self.slist.append( [net, sta, loc, chan, start, end, '', sr, '', '', '', '', '', '', '', '', ])
        #print len(self.slist)
        #exit if any channels do no contain end times
        if  null_end_time:
            sys.stderr.write("Exiting since some stations have no end time...\n")
            #should return false for automated flow system
            #caller can exit it needs
            #sys.exit(2)
            return False

    def read_stdin(self):
        #print "Using stdin"
        infh = sys.stdin
        self.populate(infh)
    def read_string(self, str):
        it = (line for line in str.split('\n'))
        self.populate(it)
    def read_string_dlt(self, str, netcode):
        it = (line for line in str.split('\n'))
        self.populate_dlt(it, netcode)
    def read_stream(self, stream):
        self.populate(stream)
    #WRITERS
    def __str__(self):
        return '\n'.join(['|'.join(line) for line in self.slist])
    def line2str(self, line):
        return  '|'.join(line) 
    def write(self, outfh):
        for line in self.slist:
            outfh.write('|'.join(line) + '\n')
    def write_file(self, filename):
        outfh = open(filename, 'w')
        self.write(outfh)
        outfh.close()
    def write_files_split(self, filename, max_lines= 500000):
        #doesn't work!! out of memory building inital sync use synclib.split_and_write
        '''With large syncs processes can run out of memory i.e. confirming
        returns list of filenames that have been written too
        '''
        filenum = 0
        start = 0
        end = max_lines
        currentfile = filename + '.%d' % filenum
        fileswritten = [currentfile]
        outfh = open(currentfile, 'w')
        while True:
            for line in self.slist[start:end]:
                outfh.write('|'.join(line) + '\n')
            laststa = self.get_sncl(line)[1]
            try:
                line = self.slist[end]
            except:
                break
            while  self.get_sncl(line)[1] == laststa:
                outfh.write('|'.join(line) + '\n')
                end += 1
                line = self.slist[end]
            outfh.close()
            start = end
            end = start + max_lines
            filenum += 1
            currentfile = filename + '.%d' % filenum
            fileswritten.append(currentfile)
            outfh = open(currentfile, 'w')
        outfh.close()
        return fileswritten
                
    def write_stdout(self):
        self.write(sys.stdout)
    def print_ms_summary(self, tolerance = 86400):
        self.sort_mash()
        if tolerance >= 129600:
            format_ts = format_day_granularity
        else:
            format_ts = format_second_granularity
        sum_dict = self.mk_sum_dict(tolerance)
        keys = sum_dict.keys()
        keys.sort()
        for key in keys:
            print "%s:%s:%s:%s:%s" % (key[0], key[1], key[2], key[3], key[4])
            for ts in sum_dict[key]:
                print format_ts(ts)
    @staticmethod
    def get_sncl(line):
        '''Note not sncl but
        net, sta, chan, loc, year, day
        for comparing to filenames
        '''
        #print self.slist
        #line = self.slist[0]
        date = line[ST].split(',')
        year = int(date[YEAR])
        day = int(date[DAY])
        return line[NET], line[STA], line[CHAN], line[LOC], year, day
    def iter_sncl(self):
        ''' returns an iterater of all sncls 1 for each sync line'''
        return itertools.imap(self.get_sncl, self.slist)

    #GENERICS
    @staticmethod
    def get_today():
        lt = localtime()
        return "%d,%d" % (lt.tm_year, lt.tm_yday)
    def keyer(self, line):
        '''returns a key for grouping synclines by sta loc chan
        STA, CHAN, LOC, SR 
        '''
        return tuple( line[:4] + [str(float(line[SR]))] )
    def keyer_by_sta(self, line):
        '''returns a key for grouping synclines by sta
        STA 
        '''
        return line[STA]
    def keyer_dlt(self, line):
        '''returns a key for grouping synclines by sta loc chan NO SR
        STA, CHAN, LOC
        '''
        return tuple( line[:4] )
    def keyer_with_day(self, line):
        '''returns a key for grouping synclines by sta loc chan day
        STA, CHAN, LOC, "(start)YEAR,DAY", SR 
        '''
        sy, sd, st = line[ST].split(',')
        return tuple( line[:4] + ["%s,%s" %(sy,sd) ,str(float(line[SR]))] )
    def mk_group(self):
        '''Returns and itertools grouby with keyer keyed'''
        return itertools.groupby( self.slist , self.keyer)
    def mk_sta_group(self):
        '''Returns and itertools grouby with keyer keyed for all channels in station'''
        return itertools.groupby( self.slist , self.keyer_by_sta)
    def mk_group_dlt(self):
        '''Returns and itertools grouby with keyer keyed'''
        return itertools.groupby( self.slist , self.keyer_dlt)
    def mk_group_with_day(self):
        '''Returns and itertools grouby with keyer keyed'''
        return itertools.groupby( self.slist , self.keyer_with_day)
    def mk_dict(self):
        '''returns a dict of {key:[syncline, syncline],..}
        in the same for as the mk_group above'''
        return dict( ( (k, list(g)) for k,g in self.mk_group()) )
    #speed testing
    def mk_epochdict_old(self):
        d = dict()
        for key, group in self.mk_group():
            if d.has_key(key):
                d[key].extend( [ self.list2epochs(line) for line in group] )
            #update
            else:
                d[key] = [ self.list2epochs(line) for line in group ]
        return d
    def mk_epochdict(self):
        d = dict()
        #make a local copy for speed up dget
        dget = d.get
        for key, group in self.mk_group():
            d[key] = dget(key, list()) + [ self.list2epochs(line) for line in group] 
        return d
    def mk_store_epochdict(self):
        ''' Creates an epochdict as above but stores it in self.epochdict for reuse.
        Created for confirming may files agains 1 dmc sync. '''
        d = dict()
        #make a local copy for speed up dget
        dget = d.get
        for key, group in self.mk_group():
            d[key] = dget(key, list()) + [ self.list2epochs(line) for line in group] 
        self.epochdict = d
    def mk_sta_epochdict_old(self):
        d = dict()
        for key, group in self.mk_sta_group():
            if d.has_key(key):
                d[key].extend( [ self.list2epochs(line) for line in group] )
            #update
            else:
                d[key] = [ self.list2epochs(line) for line in group ]
        return d
    def mk_sta_epochdict(self):
        d = dict()
        #make a local copy for speed up dget
        dget = d.get
        for key, group in self.mk_sta_group():
            d[key] = dget(key, list()) + [ self.list2epochs(line) for line in group] 
        return d
    def mk_epochdict_dlt_old(self):
        d = dict()
        for key, group in self.mk_group_dlt():
            #print key, group
            if d.has_key(key):
                d[key].extend( [ self.list2epochs(line) for line in group] )
            #update
            else:
                d[key] = [ self.list2epochs(line) for line in group ]
        return d
    def mk_epochdict_dlt(self):
        d = dict()
        #make a local copy for speed up dget
        dget = d.get
        for key, group in self.mk_group_dlt():
            d[key] = dget(key, list()) + [ self.list2epochs(line) for line in group] 
        return d
    def mk_store_epochdict_dlt(self):
        d = dict()
        #make a local copy for speed up dget
        dget = d.get
        for key, group in self.mk_group_dlt():
            d[key] = dget(key, list()) + [ self.list2epochs(line) for line in group] 
        self.epochdict_dlt = d
    def mk_sum_dict(self, tolerance, gap_dict=None, earlylate=None):
        '''returns a dict similar to mk epoch dict where timspans with a gaps less than tolerance are merged
        sorting and mashing should be done before calling!
        if gapdict is given it will be populated with gaps same format as the returned summary dict
        if earliest is given it will be updated with the earliest start time seen or earliest provided, likewise
        for latest'''
        if earlylate is None:
            earlylate = [ 0,0]
        earliest, latest = earlylate
        epochdict = self.mk_epochdict()
        if not isinstance(gap_dict, dict) :
            gap_dict = dict()
        keys = epochdict.keys()
        keys.sort()
        for key in keys:
            tslist = epochdict[key]
            gaplist = list()
            self.merge_tslist(tslist, tolerance, gaplist=gaplist)
            epochdict[key] = tslist
            gap_dict[key] = gaplist
            earliest = min(tslist[0][0], earliest or sys.maxint)
            latest = max(tslist[-1][-1], latest or 0)
        earlylate[:] = [earliest, latest]
        return epochdict
    def mk_sta_sum_dict(self, tolerance):
        '''returns a dict similar to mk epoch dict where timspans with a gaps less than tolerance are merged
        sorting and mashing should be done before calling!'''
        epochdict = self.mk_sta_epochdict()
        keys = epochdict.keys()
        keys.sort()
        for key in keys:
            tslist = epochdict[key]
            self.merge_tslist(tslist, tolerance)
            epochdict[key] = tslist
        return epochdict

    def list2epochs(self, list):
        '''from a syncline returns tuple [start, end] time as epochs'''
        return ( self.synctime2epoch(list[ST]) , self.synctime2epoch(list[ET]) )

    def synctime2epoch(self, s):
        '''returns epoch time from string s which is of the format
        YYYY,JJJ,HH:MM:SS or
        YYYY,JJJ,HH:MM:SS.ssss
        '''
        #print s
        #return int(mktime(strptime(s.split('.')[0],Sync.STP)))
        return int(timegm(strptime(s.split('.')[0],Sync.STP)))

    @staticmethod
    def epoch2string(epoch):
        FMT = '%Y:%j:%H:%M:%S'
        subsec = ('%0.4f' % (epoch % 1.0) )[1:]
        timestring = strftime(FMT,  gmtime(epoch)) + subsec
        return timestring


    def merge_tslist(self, tslist, m=1, gaplist=None):
        '''merges tslist to a tslist of overlapping time spans with margin m.
        operateds on tslist in place, no return
        if gaplist is supplied a list of tuple epochs of gaps is populated'''
        tslist.sort()
        i = 0
        while i+1 < len(tslist):
            #theres and overlap with margin
            #merge
            if tslist[i][1] + m >= tslist[i+1][0]:
                #merge the 2 with start of first(ok cuz sorted) and max of end time
                tslist[i] = ( tslist[i][0], max(tslist[i][1], tslist[i+1][1]) )
                del tslist[i+1]
            #No overlap with margin == gap record it
            else:
                if isinstance(gaplist, list):
                    gaplist.append( ( tslist[i][1], tslist[i+1][0] ) ) 
                i += 1
    def contains(self, ts, tslist):
        '''is ts in tslist?
        ts= time span; tuple of epochs
        tslist list of time spans; List of tuples of epochs
        '''
        m = self.m
        self.merge_tslist(tslist, m)
        #first sort the list
        #tslist.sort()
        ##Merge overlaping timespans with margin
        #i = 0
        #while i+1 < len(tslist):
            ##theres and overlap with margin
            #if tslist[i][1] + m >= tslist[i+1][0]:
                ##merge the 2 with start of first(ok cuz sorted) and max of end time
                #tslist[i] = ( tslist[i][0], max(tslist[i][1], tslist[i+1][1]) )
                #del tslist[i+1]
            #else:
                #i += 1
        #add margin to start and subtract from end
        tsstart = ts[0] + m
        tsend = ts[1] - m
        #compare
        for herets in tslist:
            if tsstart >= herets[0] and tsend <= herets[1]:
                #print tsstart, tsend
                #print herets
                return True
        return False
    @staticmethod
    def overlaps(ts, tslist, margin=0):
        '''returns true if ts overlaps with any in tslist'''
        tslist.sort()
        tsstart = ts[0]
        tsend = ts[1]
        for herets in tslist:
            if not ( tsend - margin <= herets[0] or  herets[1] - margin <= tsstart) :
                #overlap
                return True
        return False

    def has_duplicates(self):
        #sort first!
        #returns true if duplicate lines i.e. exact same!
        #return True in itertools.ifilter( lambda x: len(x) >1, itertools.groupby(self.slist) )
        for k,g in itertools.groupby(self.slist):
            if len( list(g) ) > 1:
                return True
        return False
    def has_1SNCLSR(self):
            # returns true if all lines are for the same SNCL-SR
        return 1 == len(self.mk_dict()) 

    #MANIPULATORS
    def remove_duplicates(self):
        self.slist =  [ grp[0] for grp in itertools.groupby(self.slist) ]
    def crush_nonwf(self):
        '''should be sorted first
        really mashes logs from syncs :: 1 line per multiple days as long as no missing days
        '''
        today = self.get_today()
        g = self.mk_group()
        self.slist = list()
        for key,group in g:
            #nonwf base on sample rate
            #convoluted for a machine with limited floating performance
            #sr =  key[kSR].split('.')[0] or 1
            #sr = int(sr)
            sr = ceil(float( key[kSR]))
            if sr == 0:
                # Non wf start crushing group
                current = group.next()
                current_start, current_end = self.list2epochs(current)
                for next in group:
                    next_start, next_end = self.list2epochs(next)
                    if next_start - current_end < LOGTOLERANCE:
                        #gap less than a day MERGE!
                        current[ET] = next[ET]
                        current[NS] = '0'
                        current[COM] = 'Modified'
                        current[DLMDCC] = today
                        #update current
                        current_end = next_end
                    else:
                        #Gap is larger than a day close current and reset
                        self.slist.append(current)
                        current = next
                        current_start, current_end = self.list2epochs(current)
                #group is empty write out last current
                self.slist.append(current)
            else:
                # sr > 0 wf don't crush just add back
                self.slist.extend(list(group))
    def mash_nonwf(self):
        '''Should be sorted first
        mashes 1ogs 1 line per day
        '''
        today = self.get_today()
        g = self.mk_group_with_day()
        self.slist = list()
        for key,group in g:
            #nonwf base on sample rate
            #convoluted for a machine with limited floating performance
            sr =  ceil(float(key[kdSR]))
            if sr == 0:
                next = group.next()
                small = next[ST]
                large = next[ET]
                for line in group:
                    if small > line[ST]: small = line[ST]
                    if large < line[ET]: large = line[ET]
                next[ST] = small
                next[ET] = large
                next[NS] = '0'
                next[COM] = "Modified"
                next[DLMDCC] = today
                self.slist.append( next )
            #wf
            else:
                self.slist.extend(list(group))
    #COMPARERS
    def has(self, os, m=1, badlist= None):
        '''os is other sync
        !!!!!RETURNS TRUE IF ALL OS IS IN SYNC!!!!!!!!!!
        returns true if ALL os is in this sync
        m = is margin of error
        badlist is a list of lines that failed
        '''
        ret = True
        self.m = m
        if not isinstance(badlist, list):
            badlist = list()
        if not isinstance(os, Sync):
            raise "%s is not a 'Sync' type instance. Is a %s" % (os, type(os))
        #make dict of self
        if not self.epochdict:
            self.mk_store_epochdict()
        heredict = self.epochdict

        #testdict = self.mk_dict()

        #iterator for other
        othergrp = os.mk_group()
        for key,group in othergrp:
            #print key
            #if it doesn't have the key::quick failure
            if not heredict.has_key(key):
                #print "<<< NOT FOUND: No day record"
                notfound = [self.line2str(g) for g in group]
                #for g in notfound: print  '\t' + g
                badlist.extend(notfound)
                #print '\n'
                ret = False
                continue
            #check each line
            hereepochlist = heredict[key]
            for line in group:
                #print line
                epochs = self.list2epochs(line)
                #print epochs
                #for i in hereepochlist:
                #    print i
                #for ps in testdict[key]:
                #    print self.line2str(ps)
                #print "\n"
                if self.contains(epochs, hereepochlist):
                    continue
                else:
                    #print "<<< NOT FOUND"
                    #print "Looking for:"
                    #print key
                    #print self.line2str(line)
                    #print "not found in:"
                    #for ps in testdict[key]:
                    #    print self.line2str(ps)
                    #print "\n"
                    badlist.append(self.line2str(line))
                    ret = False 
        return ret
    def has_overlap(self, os, margin = 0):
        '''
        returns True if self and os have any overlap
        '''
        if not self.epochdict:
            self.mk_store_epochdict()
        heredict = self.epochdict
        othergrp = os.mk_group()
        for key, group in othergrp:
            #if no sncl match quick fail
            if not heredict.has_key(key):
                #print "No SNCL", key
                continue
            hereepochlist = heredict[key]
            for line in group:
                epochs = self.list2epochs(line)
                if self.overlaps(epochs, hereepochlist, margin):
                    return True
        return False



    def has_dlt(self, os, m=1, badlist= None):
        '''os is other sync 
        this sync is made from dl_times from a dataless file
        so no sample rate comparisons
        !!!!!RETURNS TRUE IF ALL OS IS IN SYNC!!!!!!!!!!
        returns true if ALL os is in this sync
        m = is margin of error
        badlist is a list of lines that failed
        '''
        ret = True
        self.m = m
        if not isinstance(badlist, list):
            badlist = list()
        if not isinstance(os, Sync):
            raise "%s is not a 'Sync' type instance. Is a %s" % (os, type(os))
        #make dict of self
        if not self.epochdict_dlt:
            self.mk_store_epochdict_dlt()
        heredict = self.epochdict_dlt

        #testdict = self.mk_dict()

        #iterator for other
        othergrp = os.mk_group_dlt()
        for key,group in othergrp:
            #print key
            #if it doesn't have the key::quick failure
            if not heredict.has_key(key):
                #print "<<< NOT FOUND: No day record"
                notfound = [self.line2str(g) for g in group]
                #for g in notfound: print  '\t' + g
                badlist.extend(notfound)
                #print '\n'
                ret = False
                continue
            #check each line
            hereepochlist = heredict[key]
            for line in group:
                #print line
                epochs = self.list2epochs(line)
                #print epochs
                #for i in hereepochlist:
                #    print i
                #for ps in testdict[key]:
                #    print self.line2str(ps)
                #print "\n"
                if self.contains(epochs, hereepochlist):
                    continue
                else:
                    #print "<<< NOT FOUND"
                    #print "Looking for:"
                    #print key
                    #print self.line2str(line)
                    #print "not found in:"
                    #for ps in testdict[key]:
                    #    print self.line2str(ps)
                    #print "\n"
                    badlist.append(self.line2str(line))
                    ret = False 
        return ret

def c():
    a = [x for x in sys.argv]
    s1 = Sync()
    s1.read_file(a[1])
    s1.sort_mash()
    #s1.write_file(a[1] + ".msh")
    s2 = Sync()
    s2.read_file(a[2])
    s2.sort_mash()
    #s2.write_file(a[2] + '.msh')
    badlist = list()
    print s2.has(s1, m=1, badlist = badlist)
    print badlist
    #print s1.mk_dict()[('ZO', 'HVU4', '', 'BDF', '2007', '218', '100.0')]
    #print s2.mk_dict()[('ZO', 'HVU4', '', 'BDF', '2007', '218', '100.0')]

def cbl():
    a = [x for x in sys.argv]
    s2 = Sync()
    s2.read_file(a[2])
    fh = open(a[1], 'r')
    print "comparing"
    synclines = ( Sync(l) for l in fh)
    for s1 in synclines:
        if not s2.has(s1):
            print "NOT FOUND >>>"
def mash():
    s = Sync()
    s.read_file(sys.argv[1])
    s.sort_mash()
    s.write_file(sys.argv[1] + ".mshd")
def crush():
    s = Sync()
    s.read_file(sys.argv[1])
    s.sort_crush()
    s.write_file(sys.argv[1] + ".cshd")
def mash_print():
    s = Sync()
    s.read_file(sys.argv[1])
    s.sort_mash()
    print str(s)
def has_dups():
    s = Sync()
    s.read_file(sys.argv[1])
    s.sort_mash_keep_dups()
    print "File has duplicates?",  s.has_duplicates()
def dl_times2sync():
    s = Sync()
    dl_file = open(sys.argv[1])
    s.populate_dlt(dl_file, 'XA')
    print s
def dl_has_ms():
    net = 'ZW'
    badlist = list()
    dls = Sync()
    dl_file = open(sys.argv[1])
    dls.populate_dlt(dl_file, net)
    #dls.sort_mash_keep_dups()
    print dls
    mss = Sync()
    ms_file = open(sys.argv[2])
    mss.populate(ms_file)
    mss.sort_mash()
    print "dl has ms", dls.has_dlt(mss, badlist=badlist)
    print "Not found", badlist

def ds_sum():
    '''This is the function when call from the command line.
    summarizes a datasync similar to ds_times.
    -t is the tolerance of a 'continious' block of data defaults to 1 day == 86400 seconds
    '''
    import optparse
    parser = optparse.OptionParser(usage="Usage: %prog [options] datasync", version=VERSION)
    parser.add_option('-t', '--tolerance', action='store', type='int', dest='tolerance',
                      help='Maximum size (in seconds) of gaps to ignore',
                      default=86400  )
    parser.add_option('-m', '--ms', action='append',
                      help= 'MiniSEED sync file' )
    parser.add_option('-d', '--dl', action='append', 
                      help = 'DatalessSEED file. If given will be compared with MSEED sync given, and a summary of the MSEED sync will not be produced.' )
    parser.add_option('-p', '--plot', action='store_true',
                      help = 'Plot summaries' )
    parser.add_option('-s', '--station',  action="append",
                      help = 'Plot only stations matching STATION.')
    parser.add_option('-c', '--channel',  action="append",
                      help = 'Plot only channels matching CHANNEL.')
    parser.add_option('-C', '--excludechannel',  action="append",
                      help = 'Plot all channels excluding CHANNEL.')
    (options, args) = parser.parse_args()
    print options
    print args
    tolerance = options.tolerance
    if tolerance >= 129600:
        format_ts = format_day_granularity
    else:
        format_ts = format_second_granularity
    dl_sums = list()
    ms_sums = list()
    for dl in options.dl or ():
        '''Dataless sync given: compare'''
        '''tolerance command line arg is not used for dataless, always use zero and print day gran'''
        dl_sync = Sync()
        dl_sync.read_dl_file(dl)
        sum_dict = dl_sync.mk_sum_dict(0)
        dl_sums.append(sum_dict)
        keys = sum_dict.keys()
        keys.sort()
        for key in keys:
            #this works but wont work for plot needs to come out of dl_sums probably a Sync class method to filter
            if options.excludechannel:
                if key[CHAN] in options.excludechannel :
                    continue
            if options.channel:
                if key[CHAN] not in  options.channel:
                    continue
            if options.station:
                if key[STA] not in options.station:
                    continue
            print "%s:%s:%s:%s:%s" % (key[0], key[1], key[2], key[3], key[4])
            for ts in sum_dict[key]:
                print format_ts(ts)

    for ms in options.ms or ():
        '''no dataless sync given: produce summary'''
        insync = Sync()
        insync.read_file(ms)
        sum_dict = insync.mk_sum_dict(tolerance)
        ms_sums.append(sum_dict)
        keys = sum_dict.keys()
        keys.sort()
        for key in keys:
            #this works but wont work for plot needs to come out of dl_sums probably a Sync class method to filter
            if options.excludechannel:
                if key[CHAN] in options.excludechannel :
                    continue
            if options.channel:
                if key[CHAN] not in options.channel:
                    continue
            if options.station:
                if key[STA] not in options.station:
                    continue
            print "%s:%s:%s:%s:%s" % (key[0], key[1], key[2], key[3], key[4])
            for ts in sum_dict[key]:
                print format_ts(ts)


    if options.plot:
        #plot_exp_sum(ms_sums=ms_sums, dl_sums=dl_sums)
        plot_big_detailed(ms_sums=ms_sums, dl_sums=dl_sums)
        #plot_sums(ms_sums=ms_sums, dl_sums=dl_sums)

def plot_exp_sum(ms_sums=None, dl_sums= None):
    '''Plots all summaries with 1 line per station'''
    if not (ms_sums or dl_sums):
        return
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    if False: assert isinstance(plt, matplotlib.pyplot)
    fig = plt.figure()
    if False: assert isinstance(fig, matplotlib.pyplot.figure)
    ax = fig.add_subplot(111)
    #First build a comprehensive set of channel keys
    master_keys = list()
    for sums in itertools.chain(ms_sums, dl_sums):
        master_keys.extend(sums.keys())
    #uniquify the list (no dups) and sort it 
    master_keys = sorted({}.fromkeys(master_keys).keys())
    key_loc_dict = dict(zip(master_keys,range(len(master_keys))))

def plot_sums(ms_sums=None, dl_sums=None):
    '''plots  all dl_sums and ms_sums
    sums are an iterable (list or tuple) of sum_dicts'''
    if not (ms_sums or dl_sums):
        return
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    if False: assert isinstance(plt, matplotlib.pyplot)
    fig = plt.figure()
    if False: assert isinstance(fig, plt.figure)
    ax = fig.add_subplot(111)
    #First build a comprehensive set of channel keys
    master_keys = list()
    for sums in itertools.chain(ms_sums, dl_sums):
        master_keys.extend(sums.keys())
    #uniquify the list (no dups) and sort it 
    master_keys = sorted({}.fromkeys(master_keys).keys())
    key_loc_dict = dict(zip(master_keys,range(len(master_keys))))

    print "plotting dl_sums"
    for sum_dict in dl_sums:
        keys =  sum_dict.keys()
        keys.sort()
        for i, key  in enumerate(keys):
            i = key_loc_dict[key]
            for ts in sum_dict[key]:
                x = map(mdates.num2date, map(mdates.epoch2num, ts) )
                ax.plot_date(x, (i,i), xdate=True, ydate=False, lw=6, ls='-', alpha=.8, color='orange')
    #ms_sums next
    print "plotting ms_sums"
    for sum_dict in ms_sums:
        keys =  sum_dict.keys()
        keys.sort()
        for i, key  in enumerate(keys):
            i = key_loc_dict[key]
            for ts in sum_dict[key]:
                x = map(mdates.num2date, map(mdates.epoch2num, ts) )
                ax.plot_date(x, (i,i), xdate=True, ydate=False, lw=3, ls='-', alpha=.5, color='purple')
    #labels
    ax.set_yticks([float(index) for index in key_loc_dict.values()])
    ax.set_yticklabels([ ':'.join([str(x) for x in key]) for key in key_loc_dict.keys()], size='x-small')
    ax.grid(False)
    #defaultsize = fig.get_size_inches()
    #fig.set_size_inches( (defaultsize[0], len(master_keys)*1))
    plt.show()
    #max_pix = 32768
    #dpi = 100
    #num_chans = len(master_keys)
    #print num_chans
    #size = num_chans/5.0
    #print "size %(size)s" % locals()
    #fig.set_size_inches( ( size, size) )
    #fig.savefig("test.png", dpi=(dpi))
def plot_big_detailed(ms_sums=None, dl_sums=None):
    '''plots  all dl_sums and ms_sums
    sums are an iterable (list or tuple) of sum_dicts'''
    if not (ms_sums or dl_sums):
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    if False: assert isinstance(plt, matplotlib.pyplot)
    fig = plt.figure()
    if False: assert isinstance(fig, plt.figure)
    ax = fig.add_subplot(111)
    #First build a comprehensive set of channel keys
    master_keys = list()
    for sums in itertools.chain(ms_sums, dl_sums):
        master_keys.extend(sums.keys())
    #uniquify the list (no dups) and sort it 
    master_keys = sorted({}.fromkeys(master_keys).keys())
    key_loc_dict = dict(zip(master_keys,range(len(master_keys))))

    print "plotting dl_sums"
    for sum_dict in dl_sums:
        keys =  sum_dict.keys()
        keys.sort()
        for i, key  in enumerate(keys):
            i = key_loc_dict[key]
            for ts in sum_dict[key]:
                x = map(mdates.num2date, map(mdates.epoch2num, ts) )
                ax.plot_date(x, (i,i), xdate=True, ydate=False, lw=6, ls='-', alpha=.8, color='orange')
    #ms_sums next
    print "plotting ms_sums"
    for sum_dict in ms_sums:
        keys =  sum_dict.keys()
        keys.sort()
        for i, key  in enumerate(keys):
            i = key_loc_dict[key]
            for ts in sum_dict[key]:
                x = map(mdates.num2date, map(mdates.epoch2num, ts) )
                ax.plot_date(x, (i,i), xdate=True, ydate=False, lw=3, ls='-', alpha=.5, color='purple')
    #labels
    ax.set_yticks([float(index) for index in key_loc_dict.values()])
    ax.set_yticklabels([ ':'.join([str(x) for x in key]) for key in key_loc_dict.keys()], size='xx-small')
    ax.grid(False)
    #plt.show()
    #defaultsize = fig.get_size_inches()
    #fig.set_size_inches( (defaultsize[0], len(master_keys)*1))
    max_pix = 32768
    dpi = 50
    num_chans = len(master_keys)
    print num_chans
    size = num_chans/5.0
    print "size %(size)s" % locals()
    fig.set_size_inches( ( size, size) )
    fig.savefig("test.png", dpi=(dpi))
    fig.savefig("test.pdf", dpi =(dpi))

def web_plot(syncstr=None):
    '''Returns a plot for web viewing. called from view'''
    if not syncstr:
        return
    insync = Sync()
    insync.read_string(syncstr)
    ms_sums = [insync.mk_sum_dict(0),]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    if False: assert isinstance(plt, matplotlib.pyplot)
    fig = plt.figure()
    if False: assert isinstance(fig, plt.figure)
    ax = fig.add_subplot(111)
    #First build a comprehensive set of channel keys
    master_keys = list()
    for sums in itertools.chain(ms_sums):
        master_keys.extend(sums.keys())
    #uniquify the list (no dups) and sort it 
    master_keys = sorted({}.fromkeys(master_keys).keys())
    key_loc_dict = dict(zip(master_keys,range(len(master_keys))))
    #ms_sums next
    for sum_dict in ms_sums:
        keys =  sum_dict.keys()
        keys.sort()
        for i, key  in enumerate(keys):
            i = key_loc_dict[key]
            for ts in sum_dict[key]:
                x = map(mdates.num2date, map(mdates.epoch2num, ts) )
                ax.plot_date(x, (i,i), xdate=True, ydate=False, lw=3, ls='-', alpha=.5, color='purple')
    #labels
    ax.set_yticks([float(index) for index in key_loc_dict.values()])
    ax.set_yticklabels([ ':'.join([str(x) for x in key]) for key in key_loc_dict.keys()], size='small')
    ax.grid(False)
    #plt.show()
    #defaultsize = fig.get_size_inches()
    #fig.set_size_inches( (defaultsize[0], len(master_keys)*1))
    max_pix = 32768
    dpi = 60
    num_chans = len(master_keys)
    #print num_chans
    size = num_chans*2.0
    #print "size %(size)s" % locals()
    fig.set_size_inches( ( size*10, size) )
    import StringIO
    memfile = StringIO.StringIO()
    fig.savefig(memfile, dpi=(dpi))
    return memfile.getvalue()

def format_day_granularity(ts):
    '''Returns a formated string for a time span with year and julian day'''
    start = gmtime(ts[0])
    end = gmtime(ts[1])
    return '\t%0.4d:%0.3d        %0.4d:%0.3d' % ( start.tm_year, start.tm_yday, end.tm_year, end.tm_yday)
def format_second_granularity(ts):
    start = gmtime(ts[0])
    end = gmtime(ts[1])
    s1 = "%0.4d:%0.3d:%0.2d:%0.2d:%0.2d" % (start.tm_year, start.tm_yday, start.tm_hour, start.tm_min, start.tm_sec )
    s2 = "%0.4d:%0.3d:%0.2d:%0.2d:%0.2d" % (end.tm_year, end.tm_yday, end.tm_hour, end.tm_min, end.tm_sec )
    return '\t%s        %s' % (s1, s2)

def split_and_write( filename, fh, max_size = 3 * 50 * 1024**2):
    '''Takes a fh: file or stream and writes max_lines size files
    returns list of filenames'''
    assert isinstance(fh, file)
    filenum = 0
    currentfile = filename + '.%d' % filenum
    fileswritten = [currentfile]
    outfh = open(currentfile, 'w')
    while True:
        outfh.writelines(fh.readlines(max_size))
        line = fh.readline()
        #EOF
        if line == '':
            break
        sta = line.split('|')[STA]
        while sta == line.split('|')[STA] :
            outfh.write(line)
            line = fh.readline()
        if line == '':
            break
        outfh.close()
        filenum  += 1
        currentfile = filename + '.%d' % filenum
        fileswritten.append(currentfile)
        outfh = open(currentfile, 'w')
        outfh.write(line)
    outfh.close()
    return fileswritten    

### END of synclib

if __name__ == "__main__":
    ms_syncs , dl_syncs = get_command_args()
    root = Tk()
    root.geometry("%dx%d" % (WINDOW_WIDTH, WINDOW_HEIGHT) )
    root.wm_title('%s - %s' % (sys.argv[0], VERSION))
    cp = CoverPlot(root, ms_syncs, dl_syncs)
    cp.mainloop()
