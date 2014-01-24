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
	"""docstring for Reporter"""
	def __init__(self,opt,outfileDir,counter,latexTemplate="../data/template.tex",):
		self.template = open(latexTemplate,'r')
		self.docString =  TexString(string=self.template.read())
		self.outfileDir = outfileDir
		self.opt = opt
		self.counter = counter

		self.latexWriten = False
		self.imgDir = PROJECT_PATH + '/proto_err/results/%s/img/' % (self.opt.simID  )
		

	def renderOptions(self):
		outstr  = ""
		d =  vars(self.opt)
		# for k,v in d.iteritems():
		# 	tmpStr = "".join([str(k),' : ',str(v),'\\'])
		# 	outstr = outstr.join(tmpStr)
		# print outstr
		outstr = str(d)
		return outstr


	def renderTemplate(self):
		## Metadata
		self.docString.replace(k='runID',v=self.opt.simID)
		self.docString.replace(k='options',v=self.renderOptions())
		## Images
		self.docString.replace(k="SNPQualCalibration",v='%sqscoreCalibration_SNP.png' % (self.imgDir),escape=False)
		self.docString.replace(k="INSQualCalibration",v='%sqscoreCalibration_Insertion.png' % (self.imgDir),escape=False)
		self.docString.replace(k="DELQualCalibration",v='%sqscoreCalibration_Deletion.png' % (self.imgDir),escape=False)
		self.docString.replace(k="errorDistribution_hist",v='%serrorDistribution_hist.png' % (self.imgDir),escape=False)
		self.docString.replace(k="QualDistribution_hist",v='%sQualDistribution_hist.png' % (self.imgDir),escape=False)
		self.docString.replace(k="SNP_observed_simulated_expected_transition",v='%sSNP_observed_simulated_expected_transition.png' % (self.imgDir),escape=False)
		self.docString.replace(k="deletion_size_bar",v='%sdeletion_size_bar.png' % (self.imgDir),escape=False)
		self.docString.replace(k="insertion_size_bar",v='%sinsertion_size_bar.png' % (self.imgDir),escape=False)
		## Tables
		self.docString.replace(k="errorsDF",v=self.counter.errorsDF.to_latex(),escape=False)
		readCount = self.counter.readCounts.drop(['perfectAlignments','P','N','='],1)
		readCount1 = readCount[['TotalReads', 'Mapped','UnMapped','HardClippedReads']]
		readCount2 = readCount[['totalBases', 'totalAlignedBases','M','I','D','H','S']]
		self.docString.replace(k="readCount1",v=readCount1.to_latex(),escape=False)
		self.docString.replace(k="readCount2",v=readCount2.to_latex(),escape=False)

		self.docString.replace(k="significantContext",v=self.counter.contextStats[:10].to_latex(),escape=False)
		self.docString.replace(k="significantQualContext",v=self.counter.contextQualStats[:10].to_latex(),escape=False)
		self.counter.contextStats =  self.contextStats.sort(['meanQualityScoreForContext'],ascending=True)
		self.docString.replace(k="lowQualScore",v=self.counter.contextStats.to_latex(),escape=False)

	def writeLatex(self):
		self.renderTemplate()
		logging.info("Writing LaTeX to %s " % (self.outfileDir+'report.tex'))
		with open(self.outfileDir+'report.tex','w') as outfile:
			outfile.write(self.docString.text)
		self.latexWriten = True

	def generatePdfReport(self):
		if not self.latexWriten:
			self.writeLatex()
		cmd =  'pdflatex -output-directory=%s %sreport.tex  ' % (self.outfileDir,self.outfileDir)
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

		


		