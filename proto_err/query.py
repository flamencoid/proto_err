#!/usr/bin/env python
from pymongo import MongoClient
import logging
from utils import *
from pysam import AlignedRead
from error import error
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

        self.db = self.client[database]
        self.collection = collection
        self.errors = self.db[self.collection]

    def addErrors(self,errorList):
        """
        Add errors from errorList to database
        """
        self.logger.info('Uploading errors to database "proto_err.%s" for downstream querying'%(self.collection))
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
    def find(self,query,filt=None):
        if filt:
            return self.errors.find(query,filt)
        else:
            return self.errors.find(query)
    def find_errors(self,query,filt=None):
        errorList = []
        for document in self.find(query,filt):
            read = AlignedRead()
            read.seq = str(document['read'])
            errorList.append(error(true=document['true'],
                                        emission=document['emission'],
                                        read=read,readPos=document['readPos']))
        return errorList
    def deleteAll(self):
        self.logger.info("### Wiping error DB")
        self.errors.remove()
    # def addOptionsDocument(self,opt):
    #     # self.creationTime = 
    #     for key in getKeysFromValuesObject(opt):
    #         print key



