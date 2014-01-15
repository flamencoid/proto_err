## Everything to do with plotting
import math, datetime
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import logging
import json
import os
base = importr('base')
stats = importr('stats')
datasets = importr('datasets')
grdevices = importr('grDevices')
class plotter():
    """
        Plotting object from which histPlotter inherits


        Parameters
        ----------
        opt : Value
           	options passed from OptionsParser()

        """
    def __init__(self,opt):
        self.logger = logging
        self.setup(opt)
    def setup(self,opt):
       	if not opt.imgDir:
       		self.logger.error("Option must include imgDir") 
       	else:
       		self.imgDir = opt.imgDir
       	if not opt.jsonDir:
       		self.logger.error("Option must include jsonDir") 
       	else:
       		self.jsonDir = opt.jsonDir
        


class histPlotter(plotter):
    """
        Plots a histogram

        Parameters
        ----------
        dic : string
            truth base(s)
        opt : Value
           	options passed from OptionsParser()
        filename : string
            name of output

        Examples
        --------
        	>>> from optparse import Values
        	>>> opt = Values()
        	>>> opt.imgDir = '/img'
        	>>> opt.jsonDir = '/json'
            >>> plotter = histPlotter({'AAA':1,'TTT':10},opt=opt,filename="image")
            INFO:Image saved to /img/image.png
			INFO:Raw data saved to /json/image.json
        """
    def __init__(self,dic,opt,filename):
        plotter.__init__(self,opt)
        self.dic = dic
        self.imgFilename = self.imgDir + filename+ ".png"
        self.jsonFilename = self.jsonDir + filename+ ".json"
        topLevelDir = "/".join(self.imgFilename.split('/')[:-1])
        if not os.path.exists(topLevelDir):
            os.makedirs(topLevelDir)  
        topLevelDir = "/".join(self.jsonFilename.split('/')[:-1])
        if not os.path.exists(topLevelDir):
            os.makedirs(topLevelDir) 
    def plot(self):
    	"""Method to call plot"""
        dataf= ro.DataFrame({'kmer': ro.StrVector(tuple(self.dic.keys())), 
        					'Counts': ro.IntVector(tuple(self.dic.values()))} )
        gp = ggplot2.ggplot(dataf)
        pp = gp + \
             ggplot2.aes_string(x='kmer',y='Counts') + \
             ggplot2.geom_bar(stat='identity') + \
             ggplot2.theme_bw()+\
             ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90,hjust = 1)})
        grdevices.png(file=self.imgFilename, width=512, height=512)
        pp.plot()
        grdevices.dev_off()

        self.logger.info("Image saved to "+ self.imgFilename )
        with open(self.jsonFilename,'wb') as jsonOutFile:
        	json.dump(self.dic,jsonOutFile)
        self.logger.info("Raw data saved to "+ self.jsonFilename )

class densityPlotterFromLists(plotter):
    """
        Plots a histogram

        Parameters
        ----------
        dic : 
            a dictonary of lists
        opt : Value
            options passed from OptionsParser()
        filename : string
            name of output
        """
    def __init__(self,dic,opt,filename):
        plotter.__init__(self,opt)
        self.dic = dic
        self.imgFilename = self.imgDir + filename+ ".png"
        self.jsonFilename = self.jsonDir + filename+ ".json"
    def plot(self,geom='dens'):
        """Method to call plot"""
        dfDic = {}

        valueList = []
        nameList = []
        for key,values in self.dic.iteritems():
            valueList.extend(values)
            nameList.extend([key]*len(values))
        dataf= ro.DataFrame({'name': ro.StrVector(tuple(nameList)), 
                            'value': ro.FloatVector(tuple(valueList))})
        gp = ggplot2.ggplot(dataf)
        if geom == 'dens':
            pp = gp + \
                 ggplot2.aes_string(x='value',fill='factor(name)') + \
                 ggplot2.geom_density(alpha=0.5,size=2) + \
                 ggplot2.theme_bw()+\
                 ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90,hjust = 1)})
        else:
            pp = gp + \
                 ggplot2.aes_string(x='value',fill='factor(name)') + \
            ggplot2.geom_histogram(position='dodge') + \
            ggplot2.theme_bw()+\
            ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90,hjust = 1)})
        grdevices.png(file=self.imgFilename, width=512, height=512)
        pp.plot()
        grdevices.dev_off()

        self.logger.info("Image saved to "+ self.imgFilename )
        with open(self.jsonFilename+'_'+ geom,'wb') as jsonOutFile:
            json.dump(self.dic,jsonOutFile)
        self.logger.info("Raw data saved to "+ self.jsonFilename )

class multiHistPlotter(histPlotter):
    """
        Plots a histogram with multiple 

        Parameters
        ----------
        dic : string
            truth base(s)
        opt : Value
           	options passed from OptionsParser()
        filename : string
            name of output

        Examples
        --------
        	>>> from optparse import Values
        	>>> opt = Values()
        	>>> opt.imgDir = '/img'
        	>>> opt.jsonDir = '/json'
            >>> plotter = histPlotter({'AAA':{'Expected':1,'Observed':3},'TTT':{'Expected':10,'Observed':30}},,opt=opt,filename="image")
            INFO:Image saved to /img/image.png
			INFO:Raw data saved to /json/image.json
        """
    def __init__(self,dic,opt,filename):
    	histPlotter.__init__(self,dic,opt,filename)
    def plot(self):
    	"""Method to call plot"""
    	kmerList = []
    	countList = []
    	countType = []
    	for kmer,nestedDic in self.dic.iteritems():
    		for name,count in nestedDic.iteritems():
    			kmerList.append(kmer)
    			countList.append(count)
    			countType.append(name)
    	dataf= ro.DataFrame({'kmer': ro.StrVector(tuple(kmerList)), 
        					'Counts': ro.IntVector(tuple(countList)),
        					'Type': ro.StrVector(tuple(countType))} )


        gp = ggplot2.ggplot(dataf)
        pp = gp + \
             ggplot2.aes_string(x='kmer',y='Counts',fill='factor(Type)') + \
             ggplot2.geom_bar(stat='identity',position="dodge") + \
             ggplot2.theme_bw()+\
             ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 90,hjust = 1)})
        grdevices.png(file=self.imgFilename, width=512, height=512)
        pp.plot()
        grdevices.dev_off()

        self.logger.info("Image saved to "+ self.imgFilename )
        with open(self.jsonFilename,'wb') as jsonOutFile:
        	json.dump(self.dic,jsonOutFile)
        self.logger.info("Raw data saved to "+ self.jsonFilename )





        
        