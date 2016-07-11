'''
    simTester.py -- functions and classes for testing the measure and returning statistical values
'''

import numpy as np
import matplotlib.pyplot as plt
from simLearner import getData, splitData, precrec, rocCurve, simLearner

from myConstants import *
from scipy.stats import ttest_1samp, ttest_rel

  
class statTracker:
    '''
        Stores all the statistics obtained in successive tries (note: numtries = number of tries), for a given similarity measure
    '''
    def __init__(self, dB, numtries):
        self.precisions = np.zeros(numtries)                    # Values of Precision 
        self.recalls = np.zeros(numtries)                       # Values of Recall
        self.f1s = np.zeros(numtries)                           # F1 scores
        self.dB = dB                                            # Step size used to compute ROC curves
        self.areas = np.zeros(numtries)                         # Areas under roc curves
        self.numsteps = int(np.ceil(1/self.dB))                 # Number of steps used to compute ROC curves
        self.tprates = np.zeros((numtries, self.numsteps))      # True positive rates (for ROC curves)
        self.fprates = np.zeros((numtries, self.numsteps))      # False positive rates (for ROC curves)
        self.num_added = 0
        
    def addStats(self, prec, recall, rc):
        '''
            Add statistics to tracker. prec = precision, recall = recall, rc = ROC curve (true pos rates, false pos rates, area under curve)
        '''
        self.precisions[self.num_added] = prec
        self.recalls[self.num_added] = recall
        self.f1s[self.num_added] = 2 * (prec * recall)/(prec+recall)
        self.tprates[self.num_added] = rc[0]
        self.fprates[self.num_added] = rc[1]
        self.areas[self.num_added] = rc[2]
        self.num_added += 1
    
    def plot_bar(self, color, label):
        '''
            Plot ROC curve. Width = error bar size = 2*std 
        '''
        tp_means = np.zeros(self.numsteps)
        tp_stds = np.zeros(self.numsteps)
        
        for k in range(self.numsteps):
            relev_rates = self.tprates[np.logical_and(k * self.dB <= self.fprates, self.fprates < (k+1) * self.dB)]
            tp_means[k] = np.mean(relev_rates)
            tp_stds[k] = np.std(relev_rates)
        
        plt.errorbar(np.linspace(0,1,self.numsteps), tp_means, yerr=tp_stds, color=color, label = label)  

    
class simTester:
    '''
        Repeatedly tests all similarity measures over a given species' dataset
    '''
    def load(self, species, takeFull = False, C = 1, dB = 0.001):
        self.data = getData(species, takeFull=takeFull)
        self.dB = dB
        self.numfeats = self.data[0].rtable.shape[1]
        self.species = species
        self.learner = simLearner(C=C, dB = dB)
        self.learner.show_output = False
        
    def testSim(self, numtries):
        '''
            Test all similarity measures numtries times
        '''

        self.weights = np.zeros((numtries, self.numfeats))
        self.intercepts = np.zeros(numtries)
               
        self.mytr = statTracker(self.dB, numtries)
        
        self.resniktr = dict()
        self.lintr = dict()
        self.jiangtr = dict()
        self.simgictr = dict()
        
        for r in ['BIO', 'CEL', 'MOL']:
            self.resniktr[r] = statTracker(self.dB, numtries)
            self.lintr[r]=statTracker(self.dB, numtries)
            self.jiangtr[r] = statTracker(self.dB, numtries)
            self.simgictr[r] = statTracker(self.dB, numtries)
        
        for i in range(numtries):
            sdata = splitData(self.data)
            self.learner.train(sdata['trainingdata'].rtable, sdata['traininglabels'])
            stats = self.learner.test(sdata['testingdata'].rtable, sdata['testinglabels'])
            
            self.weights[i] = stats['weights']
            self.intercepts[i] = stats['intercept']
            
            self.mytr.addStats(stats['precision'], stats['recall'], stats['rc'])
            
            for r in ['BIO', 'CEL', 'MOL']:
                pr_resnik = precrec(sdata['trainingdata'].resnik[r], sdata['traininglabels'], sdata['testingdata'].resnik[r], sdata['testinglabels'])
                self.resniktr[r].addStats(pr_resnik[0], pr_resnik[1], rocCurve(sdata['testingdata'].resnik[r], sdata['testinglabels'], dB = self.dB))
            
                pr_lin = precrec(sdata['trainingdata'].lin[r], sdata['traininglabels'], sdata['testingdata'].lin[r], sdata['testinglabels'])
                self.lintr[r].addStats(pr_lin[0], pr_lin[1], rocCurve(sdata['testingdata'].lin[r], sdata['testinglabels'], dB = self.dB))
            
                pr_jiang = precrec(sdata['trainingdata'].jiang[r], sdata['traininglabels'], sdata['testingdata'].jiang[r], sdata['testinglabels'])
                self.jiangtr[r].addStats(pr_jiang[0], pr_jiang[1], rocCurve(sdata['testingdata'].jiang[r], sdata['testinglabels'], dB = self.dB))
            
                pr_simgic = precrec(sdata['trainingdata'].simgic[r], sdata['traininglabels'], sdata['testingdata'].simgic[r], sdata['testinglabels'])
                self.simgictr[r].addStats(pr_simgic[0], pr_simgic[1], rocCurve(sdata['testingdata'].simgic[r], sdata['testinglabels'], dB = self.dB))

            if (i % 100 == 0 or i < 100):
                print(str(i) + ' tries finished out of ' + str(numtries))

    def plotrcs(self, root):
        '''
            Plot all ROC curves
        '''
        self.mytr.plot_bar('b', 'Logsim')
        self.resniktr[root].plot_bar('r', 'Resnik')
        self.lintr[root].plot_bar('g', 'Lin')
        self.jiangtr[root].plot_bar('m', 'Jiang-Conrath')
        self.simgictr[root].plot_bar('y', 'SimGIC')
        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.legend(loc = 'lower right')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        
        
    def output_stats(self):
        '''
            Compute and store all statistics in a csv file. 
            
            Order: first the adaptive measure, then Resnik, Lin, Jiang, and simGIC (first all of them wrt biological process sub-ontology, 
            then all of them wrt molecular function sub-ontology, finally all of them wrt cellular component sub-ontology).
        '''
        
        stats = np.zeros((13, 8))
        pvals = np.zeros(8)
        c = 0
        
        all_learners = [self.mytr] + sum([[self.resniktr[root], self.lintr[root], self.jiangtr[root], self.simgictr[root]] for root in ['BIO', 'MOL', 'CEL']], [])
        
        best_prec = 1
        best_rec = 1
        best_f1 = 1
        best_area = 1
        
        
        
        for h in all_learners:
            stats[c,0] = np.mean(h.precisions)
            stats[c,1] = np.std(h.precisions)
            stats[c,2] = np.mean(h.recalls)
            stats[c,3] = np.std(h.recalls)
            stats[c,4] = np.mean(h.f1s)
            stats[c,5] = np.std(h.f1s)
            stats[c,6] = np.mean(h.areas)
            stats[c,7] = np.std(h.areas)
            
            if c > 0: 
                if stats[c,0] > stats[best_prec, 0]:
                    best_prec = c
                if stats[c,2] > stats[best_rec, 2]:
                    best_rec = c
                if stats[c,4] > stats[best_f1, 4]: 
                    best_f1 = c
                if stats[c, 6] > stats[best_area, 6]:
                    best_area = c
                    
            c+=1
        
        np.savetxt(statfiles[self.species], stats, delimiter='\t', header='precision \t std \t recall \t std \t F1 \t std \t ROC area \t std')
        
        pvals[0], pvals[1] = ttest_rel(self.mytr.precisions, all_learners[best_prec].precisions)
        pvals[2], pvals[3] = ttest_rel(self.mytr.recalls, all_learners[best_rec].recalls)
        pvals[4], pvals[5] = ttest_rel(self.mytr.f1s, all_learners[best_f1].f1s)
        pvals[6], pvals[7] = ttest_rel(self.mytr.areas, all_learners[best_area].areas)
        
        np.savetxt(pvalfiles[self.species], np.expand_dims(pvals, 0), delimiter='\t', header = 'precision \t p-value \t recall \t  p-value \t F1 \t p-value \t ROC area \t p-value')
        
        return stats
    
    def plot_weights(self):
        plt.plot(range(1, self.weights.shape[1]+1), np.mean(self.weights, 0))
        plt.xlim(1, self.weights.shape[1]+1)
    
def cross_train(species1, species2, numEls, numtries, dB= 0.001):
    '''
        Cross-train species1 with species2, record f1 and ROC area losses 
    '''
    data1 = getData(species1)
    data2 = getData(species2)
        
    tracker1 = statTracker(dB, numtries)
    tracker2 = statTracker(dB, numtries)
    
    learner1 = simLearner(); learner1.show_output=False
    learner2 = simLearner(); learner2.show_output=False
    
    for i in range(numtries):
        split1 = splitData(data1)
        split2 = splitData(data2)
        split3 = splitData(data2)
        
        learner1.train(split1['trainingdata'].rtable[:numEls, :], split1['traininglabels'][:numEls])
        learner2.train(split2['trainingdata'].rtable[:numEls, :], split2['traininglabels'][:numEls])
        
        stats1 = learner1.test(split2['testingdata'].rtable, split2['testinglabels'])
        tracker1.addStats(stats1['precision'], stats1['recall'], stats1['rc'])
        
        stats2 = learner2.test(split3['testingdata'].rtable, split3['testinglabels'])
        tracker2.addStats(stats2['precision'], stats2['recall'], stats2['rc'])
        
        if (i % 10 == 0 and i > 0): 
            print('.'),
    
    print 
    losses = dict()
    
    losses['f1'] = tracker1.f1s/tracker2.f1s
    losses['area'] = tracker1.areas/tracker2.areas
    
        
    return losses
        
def output_crosstrain(specieslist, numEls = 1000, numtries = 1000):
    '''Cross-train all species in the list with all species in the list. Record mean values, and p-values against null hypothesis loss=1 
        (no performance difference between training and cross-training)
    '''
    ctf1s = np.zeros((len(specieslist), len(specieslist)))
    ctf1s_pv = np.zeros((len(specieslist), len(specieslist)))
    
    ctars = np.zeros((len(specieslist), len(specieslist)))
    ctars_pv = np.zeros((len(specieslist), len(specieslist)))
    
    for s1 in range(len(specieslist)):
        for s2 in range(len(specieslist)): 
            print('Training with ' + specieslist[s1] + ', testing with ' + specieslist[s2])
            losses = cross_train(specieslist[s1], specieslist[s2], numEls, numtries)
            ctf1s[s1,s2] = np.mean(losses['f1'])
            ctf1s_pv[s1, s2] = ttest_1samp(losses['f1'], 1)[1]

            ctars[s1, s2] = np.mean(losses['area'])
            ctars_pv[s1, s2] = ttest_1samp(losses['area'], 1)[1]
    
    np.savetxt(crosstrain_f1_file, ctf1s, delimiter='\t')
    np.savetxt(crosstrain_f1_pv_file, ctf1s_pv, delimiter='\t')
    
    np.savetxt(crosstrain_roc_file, ctars, delimiter='\t')
    np.savetxt(crosstrain_roc_pv_file, ctars_pv, delimiter='\t')