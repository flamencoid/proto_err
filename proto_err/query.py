#!/usr/bin/env python
from pymongo import MongoClient
from errorCount import error
import logging
from utils import *
class errordb():
	"""docstring for db"""
	def __init__(self):
		try:
			self.client = MongoClient('localhost', 27017)
		except:
			logging.error(getDBError())
			self.client = MongoClient('localhost', 27017)
		self.db = self.client.proto_err
		self.errors = self.db.errors

			


	def addErrors(self,errorList):
		"""
		Add errors from errorList to database
		"""
		errorListOut = []
		for error in errorList:
			errorListOut.append(self.addError(error))
		return errorListOut
	def addError(self,error):
		"""
		Add an error to the database
		"""
		if not type(error) is dict:
			post = error.doc
		postID = self.errors.insert(post)
		error.dbID = postID
		return error
	def find_one(self):
		return self.errors.find_one()
	def find(self):
		return self.errors.find()
	def deleteAll(self):
		self.errors.remove()



