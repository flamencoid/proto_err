## Module for reporting error assesment
import logging
import subprocess
import shlex
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROJECT_PATH = os.path.join(BASE_DIR,os.pardir)
PROJECT_PATH = os.path.abspath(PROJECT_PATH)
print PROJECT_PATH
class Reporter(object):
	"""Generate an error report for ONT reads

	...
    Parameters
    ----------
    opt : Optparse.opt
        Options from Optparse
    outfileDir : str
    	Output file directory 
    counter : object
    	Couter object
    firstRead : object
    	arbirtrary aligned read from (s/b)amfile
    latexTemplate : str
    	Path to LaTeX report template

	"""
	def __init__(self,opt,outfileDir,imgDir,counter,basicStats=None,latexTemplate="../data/template.tex",):
		self.template = open(latexTemplate,'r')
		self.docString =  TexString(string=self.template.read())
		self.outfileDir = outfileDir
		self.opt = opt
		self.counter = counter

		self.latexWriten = False
		self.imgDir = imgDir
		self.basicStats = basicStats
		

	def renderOptions(self):
		outstr  = ""
		d =  vars(self.opt)
		# for k,v in d.iteritems():
		# 	tmpStr = "".join([str(k),' : ',str(v),'\\'])
		# 	outstr = outstr.join(tmpStr)
		# print outstr
		outstr = str(d)
		return outstr

	# def addCaptionToTable(self,tableString,captionString):
	# 	caption = "\caption{%s}"
		


	def renderTemplate(self):
		## Metadata
		self.docString.replace(k='runID',v=self.opt.runID)
		self.docString.replace(k='options',v=self.renderOptions())
		## Basic stats
		if self.basicStats:
			self.docString.replace(k='numReads',v=self.basicStats.get('readsConsidered'))
			self.docString.replace(k='totalReads',v=self.basicStats.get('totalNumberofReads'))
			self.docString.replace(k='snps',v= (float(100*self.basicStats.get('snps')) )/self.basicStats.get('readsConsidered') )
			self.docString.replace(k='ins',v=float(100*self.basicStats.get('del')) / self.basicStats.get('readsConsidered'))
			self.docString.replace(k='del',v=float(100*self.basicStats.get('ins')) / self.basicStats.get('readsConsidered'))

		## Aligned read
		self.docString.replace(k='msf',v=self.opt.outDir+"align.msf",escape=False)
		## Images
		self.docString.replace(k="SNPQualCalibration",v='%sqscoreCalibration_SNP.png' % (self.imgDir),escape=False)
		self.docString.replace(k="INSQualCalibration",v='%sqscoreCalibration_Insertion.png' % (self.imgDir),escape=False)
		self.docString.replace(k="DELQualCalibration",v='%sqscoreCalibration_Deletion.png' % (self.imgDir),escape=False)
		self.docString.replace(k="errorDistribution_hist",v='%serrorDistribution_hist.png' % (self.imgDir),escape=False)
		self.docString.replace(k="QualDistribution_hist",v='%sQualDistribution_hist.png' % (self.imgDir),escape=False)
		self.docString.replace(k="SNP_observed_expected_transition",v='%sSNP_observed_vs_expected_transition.png' % (self.imgDir),escape=False)
		self.docString.replace(k="deletion_size_bar",v='%sdeletion_size_bar.png' % (self.imgDir),escape=False)
		self.docString.replace(k="insertion_size_bar",v='%sinsertion_size_bar.png' % (self.imgDir),escape=False)
		## Tables
		self.docString.replace(k="errorsDF",v=self.counter.errorsDF.to_latex(),escape=False)
		readCount = self.counter.readCounts.drop(['perfectAlignments','P','N','='],1)
		readCount1 = readCount[['TotalReads', 'Mapped','UnMapped','HardClippedReads']]
		readCount2 = readCount[['totalBases', 'totalAlignedBases','M','I','D','H','S']]
		self.docString.replace(k="readCount1",v=readCount1.to_latex(),escape=False)
		self.docString.replace(k="readCount2",v=readCount2.to_latex(),escape=False)
		## Contexts with high significant differecne
		topTenContext = self.counter.contextStats[:10]
		topTenContext1 = topTenContext[['ContextTrue','ContextEmit','avgQual','simCount','samCount','expectedCount']]
		topTenContext2 = topTenContext[['ContextTrue','ContextEmit','pvalue','pvalue-adjust']]
		self.docString.replace(k="significantContext1",v=topTenContext1.to_latex(),escape=False)
		self.docString.replace(k="significantContext2",v=topTenContext2.to_latex(),escape=False)
		## Contexts with high significant differecne
		# topTenContext = self.counter.contextQualStats[:10]
		# topTenContext1 = topTenContext[['ContextTrue','ContextEmit','qscore','simCount','samCount','expectedCount']]
		# topTenContext2 = topTenContext[['ContextTrue','ContextEmit','qscore','pvalue','pvalue-adjust']]
		# self.docString.replace(k="significantQualContext1",v=topTenContext1.to_latex(),escape=False)
		# self.docString.replace(k="significantQualContext2",v=topTenContext2.to_latex(),escape=False)
		## Contexts with low Qscore
		topTenContext = self.counter.contextStats.sort(['avgQual'],ascending=True)[:10]
		topTenContext1 = topTenContext[['ContextTrue','ContextEmit','avgQual','simCount','samCount','expectedCount']]
		topTenContext2 = topTenContext[['ContextTrue','ContextEmit','pvalue','pvalue-adjust']]
		self.docString.replace(k="lowQualScore1",v=topTenContext1.to_latex(),escape=False)
		self.docString.replace(k="lowQualScore2",v=topTenContext2.to_latex(),escape=False)
		## Insertion Bias
		topTenContext = self.counter.contextINSStats[:10]
		topTenContext1 = topTenContext[['ContextTrue','avgQual','simCount','samCount','expectedCount','pvalue','pvalue-adjust']]
		self.docString.replace(k="inssignificantContext1",v=topTenContext1.to_latex(),escape=False)
		## Deletion Bias
		topTenContext = self.counter.contextDELStats[:10]
		topTenContext1 = topTenContext[['ContextTrue','avgQual','simCount','samCount','expectedCount','pvalue','pvalue-adjust']]
		self.docString.replace(k="delsignificantContext1",v=topTenContext1.to_latex(),escape=False)

	def writeLatex(self):
		self.renderTemplate()
		logging.info("Writing LaTeX to %s " % (self.outfileDir+'report.tex'))
		with open(self.outfileDir+'error_report.tex','w') as outfile:
			outfile.write(self.docString.text)
		self.latexWriten = True

	def generatePdfReport(self):
		if not self.latexWriten:
			self.writeLatex()
		cmd =  'pdflatex -output-directory=%s %serror_report.tex  ' % (self.outfileDir,self.outfileDir)
		print cmd
		logging.info("Running %s" % cmd)
		proc=subprocess.Popen(shlex.split(cmd))
		proc.communicate()




class TexString(object):
	"""docstring for TexString"""
	def __init__(self, string):
		self.text = str(string)
		self._latex_special_chars = {
								    '&':  r'\&',
								    '%':  r'\%',
								    '$':  r'\$',
								    '#':  r'\#',
								    '_':  r'\_',
								    '{':  r'\{',
								    '}':  r'\}',
								    '~':  r'\lettertilde{}',
								    '^':  r'\letterhat{}',
								    '\\': r'\letterbackslash{}',
								    '\n': r'\\\\',
									}

	def replace(self,k, v,escape=True):
		if escape:
			v = self.escape(v)
		self.text = self.text.replace("{{%s}}"%k,v)

	def escape(self,s):
		return ''.join(self._latex_special_chars.get(c, c) for c in s)

		


		