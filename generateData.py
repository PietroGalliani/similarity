'''
    generateData --- compute local similarity measures for training and testing logSim. Also computes non-adaptive similarity measures, for testing
'''

from myConstants import *
import MySQLdb
import numpy as np
from collections import defaultdict
import simil_data

class dataGenerator():
    '''
        This class can be used to generate (and save) data for a given species. 
    '''
    def __init__(self, species):
       
        self.species = species
        self.db = MySQLdb.connect('127.0.0.1', dbuser, dbpasswd, dbname)
        self.cursor = self.db.cursor()
        
        self.slims = self.getslims()                        #
        self.ancs, self.succs = self.get_ancs_succs()
        self.ics = self.getics()
        self.roots = self.getroots()        
    
    def savedata(self):
        '''
            Compute all the similarity values for the species (interactiong - core and full - and non-interacting pairs) and store the result
        '''
        idata = self.inter_data(False)
        idata.saveToFile(simfiles[self.species])
        idatafull = self.inter_data(True)
        idatafull.saveToFile(simfilesfull[self.species])

        nidata = self.noninter_data()
        nidata.saveToFile(nonsimfiles[self.species])
        
    def inter_data(self, useFull = False):
        '''
            Compute and return the similarity measures (including the local similarity values) for interacting pairs. If useFull, consider all 
            interacting pairs; otherwise, consider only the ones in the core dataset
        '''
        if useFull: 
            getInteracting = 'SELECT * FROM ' + self.species + '_int_full'
        else:
            getInteracting = 'SELECT * FROM ' + self.species + '_int'
        self.cursor.execute(getInteracting)
        interactingPairs = self.cursor.fetchall()
        return self.compute_simdata(interactingPairs)
    
    
    def noninter_data(self):
        '''
            Compute and return the similarity measures (including the local similarity values) for non-interacting pairs.
        '''
        getNotInteracting = 'SELECT * FROM ' + self.species+'_notint'
        self.cursor.execute(getNotInteracting)
        notInteractingPairs = self.cursor.fetchall()
        return self.compute_simdata(notInteractingPairs)
    
    
    
    def compute_simdata(self, pairsList):
        '''
            Compute all the similarities for the given list of pairs of proteins
        '''
        numPairs = len(pairsList)
        numSlims = len(self.slims)
        
        sdata = simil_data.similData()
        sdata.roots = set(self.roots)
        sdata.num_elements = numPairs
        sdata.rtable = np.zeros((numPairs, numSlims))
        sdata.relics = np.zeros(numPairs)
        
        for r in self.roots: 
            sdata.resnik[r] = np.zeros(numPairs)
            sdata.lin[r] = np.zeros(numPairs)
            sdata.jiang[r] = np.zeros(numPairs)
            sdata.simgic[r] = np.zeros(numPairs)
            
        
        c = 0
        for p in pairsList:
            a = self.get_annots(p[0])
            a_dir = self.get_direct_annots(p[0])
            b = self.get_annots(p[1])   
            b_dir = self.get_direct_annots(p[1])
            
            sdata.rtable[c, :] = self.resnik_slim_aux(a, b)
            
            sdata.relics[c] = self.relics_aux(a_dir, b_dir)
            
            [m_ic, p1_ic, p2_ic] = self.MICA_aux(a_dir, b_dir)
            simgic = self.simGIC_aux(a, b)

            for r in self.roots: 
                sdata.resnik[r][c] = sdata.rtable[c, self.roots[r][0]]
                
                if (m_ic[r] > 0):
                    sdata.lin[r][c] = 2* m_ic[r]/(p1_ic[r] + p2_ic[r])
                sdata.jiang[r][c] = 1+ m_ic[r] - (p1_ic[r] + p2_ic[r])/2#1/(1 + p1_ic[r] + p2_ic[r] - 2 * m_ic[r])
                sdata.simgic[r][c] = simgic[r]

            c+=1
        return sdata
        
    def get_annots(self, prod):
        '''
            Get all the GO terms that are associated -- directly, or indirectly throught their successors -- to the given product 
        '''
        getAncestorsCommand = 'SELECT DISTINCT anc.id\
            FROM term as anc INNER JOIN graph_path ON anc.id = graph_path.term1_id\
            INNER JOIN term on term.id = graph_path.term2_id\
            INNER JOIN '+self.species+'_assoc on ' + self.species+'_assoc.go_id = term.id\
            WHERE prod = %s'
        
        self.cursor.execute(getAncestorsCommand, {prod})
        annots = self.cursor.fetchall()
        return set([a[0] for a in annots])
    
    def get_direct_annots(self, prod):
        '''
            Get all the GO terms that are _directly_ associated to the given product
        '''
        getAncestorsCommand = 'SELECT DISTINCT term.id FROM term\
            INNER JOIN ' + self.species + '_assoc on ' + self.species+'_assoc.go_id = term.id WHERE prod = %s'
            
        
        self.cursor.execute(getAncestorsCommand, {prod})
        annots = self.cursor.fetchall()
        return set([a[0] for a in annots])
    
    def get_ancs_succs(self):
        '''
            For any GO term, get the list of all its ancestors and successors
        '''
        
        getPaths = 'SELECT term1_id, term2_id FROM graph_path AS gp INNER JOIN term AS rel ON rel.id = gp.relationship_type_id'# WHERE (rel.name = \'part_of\' OR rel.name = \'is_a\')'
        self.cursor.execute(getPaths)
        myPaths = self.cursor.fetchall()
        
        succDict = defaultdict(set)
        ancDict = defaultdict(set)
        for p in myPaths:
            succDict[p[0]].add(p[1])
            ancDict[p[1]].add(p[0])
        
        
        return (ancDict, succDict)
    
    def getslims(self):
        '''
            Get the slim sub-ontology
        '''
        self.cursor.execute('select slim_id, go_id from generic_slim order by slim_id')
        slims = dict(self.cursor.fetchall())
        return slims
    
    def getics(self):
        '''
            Get the information contents
        '''
        self.cursor.execute('select term_id, ic from ' + self.species+'_ic')
        return defaultdict(int, self.cursor.fetchall())
            
    def getroots(self):
        '''
            For any root of the Gene Ontology ('BIO' = biological process, 'MOL' = molecular function, 'CEL' = cellular component, 
            store the corresponding index of every term in the slim array and the corresponding GO id number 
        '''
        roots = dict()
        getslimid = 'SELECT slim_id - 1, go_id FROM generic_slim INNER JOIN term ON term.id = go_id WHERE name = %s'
        self.cursor.execute(getslimid, {'biological_process'})
        roots['BIO'] = self.cursor.fetchall()[0]
        self.cursor.execute(getslimid, {'molecular_function'})
        roots['MOL'] = self.cursor.fetchall()[0]
        self.cursor.execute(getslimid, {'cellular_component'})
        roots['CEL'] = self.cursor.fetchall()[0]
        return roots
    
    
    def resnik_slim_aux(self, annots1, annots2):
        '''
            Compute the partial similarity measures between two (upwards-closed) sets of annotations
        '''
        common_annots = annots1.intersection(annots2)
        
        res_ics = np.zeros(len(self.slims))
        
        
        for s in self.slims: 
            shared_sub_s = common_annots.intersection(self.succs[self.slims[s]])
            if len(shared_sub_s) > 0: 
                res_ics[s-1] = np.max([self.ics[i] for i in shared_sub_s])
           
        return res_ics
    
    def relics_aux(self, a1, a2):
        '''
            Compute the total information content of the two sets of annotations
        '''
        return sum([self.ics[a] for a in a1 if a in self.ics]) + sum ([self.ics[a] for a in a2 if a in self.ics]+[0])
      
    def MICA_aux(self, a1, a2):
        '''
            Given two sets of annotations, corresponding to the terms _directly_ associated to two proteins, compute the IC of the most informative 
            common ancestors in any branch of the gene ontology. Return also the ICs of the corresponding terms in the two annotation sets. 
        '''
        m_ic = defaultdict(int)
        p1_ic = defaultdict(int)
        p2_ic = defaultdict(int)
        
        #t1_name = ''
        #t2_name = ''
        #mica_name = ''
        
        for r in self.roots:
            for a in a1.intersection(self.succs[self.roots[r][1]]):
                for b in a2.intersection(self.succs[self.roots[r][1]]):
                    common_ancs = self.ancs[a].intersection(self.ancs[b])
                    mica_ic = max([self.ics[c] for c in common_ancs])
                    if mica_ic > m_ic[r]: 
                        m_ic[r] = mica_ic
                        p1_ic[r] = self.ics[a]
                        p2_ic[r] = self.ics[b]
                        #t1_name = a
                        #t2_name = b
                        #mica_name = [c for c in common_ancs if self.ics[c] == mica_ic][0]
        
            #print(str(t1_name) + ', ' + str(t2_name) + ': ' + str(mica_name))           
        return [m_ic, p1_ic, p2_ic]
    
    def simGIC_aux(self, a1, a2):
        '''
            Given two upwards-closed sets of annotations, corresponding to the GO terms associated - even indirectly - to two proteins, 
            compute the value of the simGIC similarity wrt every root of the ontology.
        '''
        i = a1.intersection(a2)
        u = a1.union(a2)
        simgic = dict()
        for r in self.roots: 
            ir = i.intersection(self.succs[self.roots[r][1]])
            ur = u.intersection(self.succs[self.roots[r][1]])
            numerator = sum([self.ics[t] for t in ir])
            if (numerator == 0): 
                simgic[r] = 0
            else:
                denominator = sum([self.ics[t] for t in ur])
                simgic[r] = numerator/denominator
        return simgic
    
def generateAll():
    '''
        Computes and stores all similarity data for all species
    '''
    for species in allspecies: 
        print(species + ': computing similarity measures...')
        d = dataGenerator(species)
        d.savedata()
