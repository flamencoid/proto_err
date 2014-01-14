#!/usr/bin/env python
from pymongo import MongoClient
import logging
from utils import *
from pysam import AlignedRead
from error import error
import datetime
class errordb():
    """docstring for db"""
    def __init__(self,database='proto_err',collection='errors'):
        self.logger = logging
        try:
            self.client = MongoClient('localhost', 27017)
            self.logger.info("Connection to localhost DB succesful")
        except:
            try:
                self.client = MongoClient('88.80.186.60', 27017)
                self.logger.info("Connection to DB at 88.80.186.60 succesful")
            except:
                MongoClient('localhost', 27017)
                logging.error(getDBError())
        self.database = database
        self.db = self.client[database]
        self.collection = collection
        self.errors = self.db[self.collection]
        

    def addErrors(self,errorList):
        """
        Add errors from errorList to database
        """

        self.logger.info('Uploading errors to database "%s.%s" for downstream querying'%(self.database,self.collection))
        self.errors = self.db[self.collection]
        posts = [error.doc for error in errorList]
        if posts:
            self.errors.insert(posts)

    def addError(self,error):
        """
        Add an error to the database
        """
        self.errors = self.db[self.collection]
        if not type(error) is dict:
            post = error.doc
        postID = self.errors.insert(post)
        error.dbID = postID
        return error
    def find_one(self,query,filt=None):
        self.errors = self.db[self.collection]
        if filt:
            return self.errors.find_one(query,filt)
        else:
            return self.errors.find_one(query)
    def find(self,query,filt=None):
        self.errors = self.db[self.collection]
        if filt:
            return self.errors.find(query,filt)
        else:
            return self.errors.find(query)
    def find_errors(self,query,filt=None):
        errorList = []
        for document in self.find(query,filt):
            read = AlignedRead()
            # read.seq = str(document['read'])
            read.seq = ''
            errorList.append(error(true=document['true'],
                                        emission=document['emission'],
                                        read=read,readPos=document['readPos']))
        return errorList
    def deleteAll(self):
        self.logger.info("### Wiping error DB")
        self.errors.remove()

    def addMetaData(self,opt,t,errorBias=None):
        self.md = self.db['metaData']
        self.logger.info('Uploading %s metaData to database "%s.metaData"'%(t,self.database))
        document = vars(opt)
        document['type'] = t
        document['date'] = str(datetime.date.today())
        document['time'] = str(datetime.datetime.today().hour) + ':'+str(datetime.datetime.today().minute)
        if self.md.find({'type':t}).count() > 0:
            self.logger.warning("### Metadata already found removing and re uploading")
            self.md.remove({'type':t})
        self.md.insert(document) 
        
        if errorBias:
            errorBiasDocument = {}
            for pattern,tup in errorBias.iteritems():
                errorBiasDocument[pattern.pattern] = tup
            errorBiasDocument['type'] = 'errorBias'
            if self.md.find({'type':'errorBias'}).count() > 0:
                self.logger.warning("### Errorbias metadata already found removing and re uploading")
                self.md.remove({'type':'errorBias'})
            self.md.insert(errorBiasDocument)
    def getMetaData(self):
        """
        Return the metadata associated with current database
        """
        results = []
        for res in self.db['metaData'].find():
            results.append(res)
        return results






