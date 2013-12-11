## Everything to do with plotting
import math, datetime
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import logging
base = importr('base')
stats = importr('stats')
datasets = importr('datasets')
grdevices = importr('grDevices')
class plotter():
    """plotting object"""
    def __init__(self,opt):
        self.outdir = opt.imgDir   
        self.logger = logging   
        


class histPlotter(plotter):
    """docstring for histPlotter"""
    def __init__(self,dic,opt,filename):
        plotter.__init__(self,opt)
        self.__dic = dic
        self.filename = self.outdir + filename+ ".png"
    def plot(self):
        dataf= ro.DataFrame({'kmer': ro.StrVector(tuple(self.__dic.keys())), 
        					'Counts': ro.IntVector(tuple(self.__dic.values()))} )
        gp = ggplot2.ggplot(dataf)
        pp = gp + \
             ggplot2.aes_string(x='kmer',y='Counts') + \
             ggplot2.geom_bar()
        grdevices.png(file=self.filename, width=512, height=512)
        pp.plot()
        grdevices.dev_off()
        self.logger.info("File saved to "+ self.filename )



        
        