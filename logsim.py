'''
logsim.py --- some high-level functions for loading all the data and testing the similarity measure
'''

import myConstants
import preprocess
import generateData
import simTester
import csv

def prepareAll():
    '''
        Does all the preprocessing (including generating and storing the local similarity values)
    '''
    print('Loading ontology and species-specific files...')
    preprocess.loadall()
    print('Generating partial and non-adaptive similarities...')
    generateData.generateAll()
    print('Done!')
    
def testSpecies(species, numtries = 1000):
    '''
        Tests the similarity measure over the given species. Stores the statistics in the corresponding files, and returns the myTester instance
        in case one wants to make other analyses
    '''
    if species not in myConstants.allspecies:
        print('Species not recognized. Known species: ' + str(myConstants.allspecies))
    myTester = simTester.simTester()
    print('Loading species data...')
    myTester.load(species)
    print('Testing similarity measure...')
    myTester.testSim(numtries)
    print('Saving data...')
    myTester.output_stats()
    print('Done. Means and stds are in ' + myConstants.statfiles[species] + ', p-values are in ' + myConstants.pvalfiles[species])
    
    return myTester 

def crossTrain(numEls = 1000, numtries = 1000):
    '''
        Cross-trains all species against all species, saves in the corresponding files. 
    '''
    simTester.output_crosstrain(['scere', 'hsap', 'ecoli', 'mus'], numEls=numEls, numtries=numtries)
    
def tolatex(species):
    myfile = open(myConstants.statfiles[species], 'r')
    myreader = csv.reader(myfile, delimiter='\t')
    rownames = ['logSim', 'Resnik (BIO)', 'Lin (BIO)', 'Jiang (BIO)', 'simGIC (BIO)', 
                          'Resnik (MOL)', 'Lin (MOL)', 'Jiang (MOL)', 'simGIC (MOL)',
                          'Resnik (CEL)', 'Lin (CEL)', 'Jiang (CEL)', 'simGIC (CEL)']
    
    c = 0
    myreader.next()
    for l in myreader:
        print(rownames[c] + ' & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f &  %.4f $\pm$ %.4f & %.4f $\pm$ %.4f\\\\' % tuple([float(v) for v in l]))
        
        c+=1
    