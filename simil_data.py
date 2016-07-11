'''
    simil_data.py -- a simple class for storing and manipulating the similarity data 
'''
from myConstants import *
import numpy as np
import pickle


class similData:
    def __init__(self):
        self.roots = None           # The roots of the ontology
        self.rtable = None          # The local similarities
        self.relics = None          # The total relation information contents for the pairs
        self.num_elements = 0       # The number of pairs in the data

                                    # For all roots (and all pairs):
        self.resnik = dict()        # The Resnik similarities
        self.lin = dict()           # The Lin similarities
        self.jiang = dict()         # The Jiang-Conrath similarities
        self.simgic = dict()        # The SimGIC similarities
        
    def saveToFile(self, filename):
        '''
            Store the similarity data in the file
        '''
        with open(filename, 'wb') as outp: 
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)
            
    @staticmethod
    def loadFromFile(filename):
        '''
            Load the similarity data from a file
        '''
        with open(filename, 'rb') as inp: 
            return pickle.load(inp)
        
    def shuffle(self, svec):
        '''
            Rearrange all the data with respect to the given permuation
        '''
        self.rtable = self.rtable[svec, :]
        self.relics = self.relics[svec]
        for r in self.roots:
            self.resnik[r] = self.resnik[r][svec]
            self.lin[r] = self.lin[r][svec]
            self.jiang[r] = self.jiang[r][svec]
            self.simgic[r] = self.simgic[r][svec]
            
    def select(self, svec):
        '''
            Select a subset of the data. Svec lists the indexes of the pairs to include
        '''
        newdata = similData()
        newdata.roots = self.roots
        newdata.num_elements = len(svec)
        newdata.rtable = self.rtable[svec, :]
        newdata.relics = self.relics[svec]
        for r in self.roots: 
            newdata.resnik[r] = self.resnik[r][svec]
            newdata.lin[r] = self.lin[r][svec]
            newdata.jiang[r] = self.jiang[r][svec]
            newdata.simgic[r] = self.simgic[r][svec]
    
        return newdata
    
    def concatenate(self, otherData):
        '''
            Concatenate this similarity dataset with another one. Return the result. 
        '''
        newdata = similData()
        newdata.num_elements = self.num_elements + otherData.num_elements
        newdata.roots = self.roots
        newdata.rtable = np.concatenate((self.rtable, otherData.rtable), 0)

        newdata.relics = np.concatenate((self.relics, otherData.relics))
        
        for r in self.roots:
            newdata.resnik[r] = np.concatenate((self.resnik[r], otherData.resnik[r]))
            newdata.lin[r] = np.concatenate((self.lin[r], otherData.lin[r]))
            newdata.jiang[r] = np.concatenate((self.jiang[r], otherData.jiang[r]))
            newdata.simgic[r] = np.concatenate((self.simgic[r], otherData.simgic[r]))
        
        return newdata        