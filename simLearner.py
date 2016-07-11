'''
    simLearner.py -- functions and classes for training and evaluating the logSim measure (and computing the others, for comparison) 
'''

from sklearn.linear_model import LogisticRegression
import numpy as np
import matplotlib.pyplot as plt
from simil_data import similData
from myConstants import *
from scipy import stats

  
def prepareData(interactingPairs, notInteractingPairs, interactingRelICs, notInteractingRelICs):    
    '''
        Combine interacting and non-interacting data into a single dataset. Return the data, the labels, and the relation ICs 
    '''
    numIntPairs = interactingPairs.shape[0]
    numNotIntPairs = notInteractingPairs.shape[0]
        
    allData = np.concatenate((interactingPairs, notInteractingPairs), 0)
    allLabels = np.zeros(numIntPairs + numNotIntPairs); allLabels[:numIntPairs] = 1
    allRelICs = np.concatenate((interactingRelICs, notInteractingRelICs), 0)
    
    dataShuffle = np.arange(numIntPairs + numNotIntPairs); np.random.shuffle(dataShuffle)
    allData = allData[dataShuffle, :]
    allLabels = allLabels[dataShuffle]
    allRelICs = allRelICs[dataShuffle]
    
    data = dict()
    data['data'] = allData
    data['labels'] = allLabels
    data['relIC'] = allRelICs
    return data

def splitData(allData):
    '''
        Split the given data into a training and a testing set.
    '''    
    idata = allData[0]
    numIntPairs = idata.num_elements
    
    nidata = allData[1]
    numNotIntPairs = nidata.num_elements
    
    numpergroup = min(numIntPairs, numNotIntPairs)
    if (numpergroup % 2 != 0):
        numpergroup -= 1
    trainingSize = numpergroup
    testingSize = numpergroup
    
    intPairs_training = np.random.choice(numIntPairs, trainingSize/2, replace=False)
    nonIntPairs_training = np.random.choice(numNotIntPairs, trainingSize/2, replace=False)
    
    available_int_testing = np.delete(np.arange(numIntPairs), intPairs_training)
    available_notint_testing = np.delete(np.arange(numNotIntPairs), nonIntPairs_training)
    
    intPairs_testing = np.random.choice(available_int_testing, numpergroup/2, replace=False)
    nonIntPairs_testing = np.random.choice(available_notint_testing, numpergroup/2, replace=False)
    
    training_idata = idata.select(intPairs_training)
    training_nidata = nidata.select(nonIntPairs_training)
    trainingData = training_idata.concatenate(training_nidata)
    trainingLabels = np.zeros(trainingSize); trainingLabels[:trainingSize/2] = 1;
    
    trainingShuffle = np.arange(trainingSize); np.random.shuffle(trainingShuffle); 
    
    trainingData.shuffle(trainingShuffle)
    trainingLabels = trainingLabels[trainingShuffle]; 
    
    
    
    testingData = idata.select(intPairs_testing).concatenate(nidata.select(nonIntPairs_testing))
    testingLabels = np.zeros(testingSize); testingLabels[:testingSize/2] = 1;
    
    testingShuffle = np.arange(testingSize); np.random.shuffle(testingShuffle); 
    
    testingData.shuffle(testingShuffle)
    testingLabels = testingLabels[testingShuffle]; 
    
    
    splitdata = dict()
    splitdata['trainingdata'] = trainingData
    splitdata['traininglabels'] = trainingLabels
    splitdata['testingdata'] = testingData
    splitdata['testinglabels'] = testingLabels
    
    return splitdata
    
def getData(species, takeFull=False):
    '''
        Load the similarity data from file. If takeFull=True, consider all interacting pairs; otherwise, consider only the core (trusted) ones. 
    '''
    if takeFull: 
        idata = similData.loadFromFile(simfilesfull[species])
    else:
        idata = similData.loadFromFile(simfiles[species])
    nidata = similData.loadFromFile(nonsimfiles[species])
    
    return (idata, nidata)

def rocCurve(sims, labels, dB = 0.00001):
    '''
        Takes a list of similarity and the corresponding labels (1 = interacting, 0 = non-interacting).
        Computes the area under the ROC curve, as well as the curve itself
    '''
    
    boundaries = np.linspace(0, 1, np.round(1/dB))
    truepos = np.zeros(len(boundaries))
    falsepos = np.zeros(len(boundaries))
    
    numpos = np.sum(labels)
    numneg = len(labels) - numpos
    area = 0
    
    for c in range(len(boundaries)):
        truepos[c] = np.sum(labels[sims > boundaries[c]])/(numpos+0.0); 
        falsepos[c] = np.sum((1-labels)[sims > boundaries[c]])/(numneg + 0.0)
        if c > 0: 
            area += (falsepos[c-1] - falsepos[c])* (truepos[c] + truepos[c-1])/2
        
    return (truepos, falsepos, area)


def precrec(trainingSims, trainingLabels, testingSims, testingLabels):
    '''
        Compute the precision and the recall for a given (non-adaptive) similarity measure. First use the training similarities and labels 
        to select - through logistic regression - the boundary between interacting and non-interacting pairs, then test this value wrt 
        the similarity values and the labels of the testing set. 
    '''
    l = LogisticRegression()
    l.fit(np.expand_dims(trainingSims,1), trainingLabels)
    predictions = l.predict(np.expand_dims(testingSims,1))
    truepos = np.sum(predictions[testingLabels==1])
    falsepos = np.sum(predictions[testingLabels == 0])
    falseneg = np.sum(1-predictions[testingLabels == 1])
        
    precision = (truepos + 0.0)/(truepos + falsepos)
    recall = (truepos + 0.0)/(truepos + falseneg) 
    
    return (precision, recall)

class simLearner:
    '''
        The adaptive similarity measure
    '''
    def __init__(self, title='Interacting and non-interacting pairs', C = 1, dB = 0.001):
        self.learner = LogisticRegression(C=C)#SVC(kernel='linear', probability=True, C=1.5)
        self.xlabel='similarity'
        self.ylabel='# occurrences'
        self.title= title
        self.show_output = True
        self.dB = dB
                
    def train(self, trainingData, trainingLabels):
        '''
            Fit the measure to the training data and lavels
        '''
        self.learner.fit(trainingData, trainingLabels)
        

    def sim(self, simData):
        '''
            Compute the similarity values (that is, the probabilities of interactiong) for the local similarity data corresponding to a list of value
        '''
        proba = self.learner.predict_proba(simData)
        return proba[:, 1]
    
    
    def cmpweights(self):
        '''
            Return the weights (for the slim ids) estimated during the training phase
        '''
        return self.learner.coef_
    
    def cmpint(self):
        '''
            Return the intercept of the underlying logistic regression model
        '''
        return self.learner.intercept_
    
    def test(self, testingData, testingLabels):
        '''
            Test the performance of the measure over a testing dataset and the corresponding labels. If show_output = true, plot the resulting similarity 
            histogram for the interacting and non-interacting pairs; otherwise, return the resulting data (weights, intercept, performance measures) 
            inside a python dictionary
        '''
        probs = self.sim(testingData)
        predictions = probs > 0.5
        truepos = np.sum(predictions[testingLabels==1])
        falsepos = np.sum(predictions[testingLabels == 0])
        falseneg = np.sum(1-predictions[testingLabels == 1])
        
        precision = (truepos + 0.0)/(truepos + falsepos)
        recall = (truepos + 0.0)/(truepos + falseneg)
        if self.show_output:
            print('precision: ' + str(precision))
            print('recall: ' + str(recall))
        
            plt.figure(1)
            plt.hist([probs[testingLabels==1], probs[testingLabels==0]], color=['green', 'red'], label=['Interacting', 'Non-interacting'] )
            plt.title(self.title)
            plt.xlabel('Similarity')
            plt.legend()
            plt.show()
            
        else: 
            stats = dict()
            stats['rc'] = rocCurve(probs, testingLabels, self.dB)
            stats['precision'] = precision
            stats['recall'] = recall
            stats['weights'] = self.cmpweights()
            stats['intercept'] = self.cmpint()
            return stats
        


class statTracker:
        def __init__(self, dB, numtries):
            self.precisions = np.zeros(numtries)
            self.recalls = np.zeros(numtries)
            self.f1s = np.zeros(numtries)
            self.dB = dB
            self.areas = np.zeros(numtries)
            self.numsteps = int(np.ceil(1/self.dB))
            self.tprates = np.zeros((numtries, self.numsteps))
            self.fprates = np.zeros((numtries, self.numsteps))
            self.num_added = 0
            
        def addStats(self, prec, recall, rc):
            self.precisions[self.num_added] = prec
            self.recalls[self.num_added] = recall
            self.f1s[self.num_added] = 2 * (prec * recall)/(prec+recall)
            self.tprates[self.num_added] = rc[0]
            self.fprates[self.num_added] = rc[1]
            self.areas[self.num_added] = rc[2]
            self.num_added += 1
        
        def plot_bar(self, color, label):
            tp_means = np.zeros(self.numsteps)
            tp_stds = np.zeros(self.numsteps)
            
            for k in range(self.numsteps):
                relev_rates = self.tprates[np.logical_and(k * self.dB <= self.fprates, self.fprates < (k+1) * self.dB)]
                tp_means[k] = np.mean(relev_rates)
                tp_stds[k] = np.std(relev_rates)
            
            plt.errorbar(np.linspace(0,1,self.numsteps), tp_means, yerr=tp_stds, color=color, label = label)  
            
        def get_precision(self):
            return (np.mean(self.precisions), np.std(self.precisions))
    
        def get_recall(self):
            return (np.mean(self.recalls), np.std(self.recalls))      

        
    
    
class simTester:
    def load(self, species, takeFull = False, C = 1, dB = 0.001):
        self.data = getData(species, takeFull=takeFull)
        self.dB = dB
        self.numfeats = self.data[0].rtable.shape[1]
        self.species = species
        self.learner = simLearner(C=C, dB = dB)
        self.learner.show_output = False
        
    def testSim(self, numtries):
        #self.precisions = np.zeros(numtries)
        #self.recalls = np.zeros(numtries)
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
            #self.precrec.addprecrec(stats['precision'], stats['recall'])
            #self.precisions[i] = stats['precision']
            #self.recalls[i] = stats['recall']
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
            
            
            
    def plotboxweights(self, title):
        fig = plt.figure()
        plt.boxplot(self.weights)
        frame = plt.gca()
        frame.axes.get_xaxis().set_ticks([])
        plt.title(title)
                
    def plotrcs(self, root):
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
        
        
    def plot_stats(self, title):
        plt.figure()
        plt.hist(self.precrec.precisions)
        plt.title(title + ', precision: mean = ' + str(np.mean(self.precrec.precisions))+ ', std = ' + str(np.std(self.precrec.precisions)))
        plt.figure()
        plt.hist(self.recalls)
        plt.title(title + ', recall: mean = ' + str(np.mean(self.precrec.recalls))+ ', std = ' + str(np.std(self.precrec.recalls)))
        plt.figure()
        plt.hist(self.intercepts)
        plt.title(title + ', intercept: mean = ' + str(np.mean(self.intercepts))+ ', std = ' + str(np.std(self.intercepts)))
    
    def plot_weight(self, title, weightid):
        plt.figure()
        plt.hist(self.weights[:, weightid])
        plt.title(title + ', weight ' + str(weightid) + ': mean = ' + str(np.mean(self.weights[:, weightid])) + ', recall = ' + np.mean(self.weights[:, weightid]))
        
    def output_stats(self):
        stats = np.zeros((13, 8))
        c = 0
        for h in [self.mytr] + sum([[self.resniktr[root], self.lintr[root], self.jiangtr[root], self.simgictr[root]] for root in ['BIO', 'MOL', 'CEL']], []):
            stats[c,0] = h.get_precision()[0]
            stats[c,1] = h.get_precision()[1]
            stats[c,2] = h.get_recall()[0]
            stats[c,3] = h.get_recall()[1]
            stats[c,4] = np.mean(h.areas)
            stats[c,5] = np.std(h.areas)
            f1scores = 2*(h.precisions*h.recalls)/(h.precisions+h.recalls)
            stats[c,6] = np.mean(f1scores)
            stats[c,7] = np.std(f1scores)
            c+=1
        np.savetxt('stats.csv', stats, delimiter=',')
        return stats
    
    def learningCurves(self, stepsize = 100):
        sdata = splitData(self.data)
        sl = simLearner()
        sl.show_output = False
        num_steps = len(range(stepsize, sdata['trainingdata'].num_elements, stepsize))
        
        precisions = np.zeros((num_steps, 2))
        recalls = np.zeros((num_steps, 2))
        areas = np.zeros((num_steps, 2))
        
        c = 0
        for i in range (stepsize, sdata['trainingdata'].num_elements, stepsize): 
            sl.train(sdata['trainingdata'].rtable[:i, :], sdata['traininglabels'][:i])
            stats_training = sl.test(sdata['trainingdata'].rtable[:i, :], sdata['traininglabels'][:i])
            stats_testing = sl.test(sdata['testingdata'].rtable, sdata['testinglabels'])
            
            precisions[c, 0] = stats_training['precision']
            recalls[c, 0] = stats_training['recall']
            areas[c, 0] = stats_training['rc'][2]
            
            precisions[c, 1] = stats_testing['precision']
            recalls[c, 1] = stats_testing['recall']
            areas[c, 1] = stats_testing['rc'][2]
            
            c+=1
            
        return (range(stepsize, sdata['trainingdata'].num_elements, stepsize), precisions, recalls, areas)
        
    def plotLearningCurves(self, stepsize=100):
        numtraining, precisions, recalls, areas = self.learningCurves(stepsize)
        plt.figure()
        plt.subplot(1,3,1)
        plt.plot(numtraining, precisions, label=['training', 'testing'])
        plt.xlabel('training set size')
        plt.ylabel('precision')
        plt.subplot(1,3,2)
        plt.plot(numtraining, recalls, label=['training', 'testing'])
        plt.xlabel('training set size')
        plt.ylabel('recall')
        plt.subplot(1,3,3)
        plt.plot(numtraining, areas, label=['training', 'testing'])
        plt.xlabel('training set size')
        plt.ylabel('ROC area')
        
def computef1s(st):
    f1s = np.zeros((13, 2))
    c = 0
    for h in [st.mytr] + sum([[st.resniktr[root], st.lintr[root], st.jiangtr[root], st.simgictr[root]] for root in ['BIO', 'MOL', 'CEL']], []):
        f1scores = 2*(h.precisions*h.recalls)/(h.precisions+h.recalls)
        f1s[c,0] = np.mean(f1scores)
        f1s[c,1] = np.std(f1scores)
        c+=1
    np.savetxt('f1s.csv', f1s, delimiter=',')
    return f1s

def f1scores(tr):
    return 2*(tr.precisions*tr.recalls)/(tr.precisions+tr.recalls)

def cross_train(species1, species2, numEls, numtries, dB= 0.001):
    data1 = getData(species1)
    data2 = getData(species2)
        
    tracker1 = statTracker(dB, numtries)
    tracker2 = statTracker(dB, numtries)
    
    learner1 = simLearner(); learner1.show_output=False
    learner2 = simLearner(); learner2.show_output=False
    
    for _ in range(numtries):
        split1 = splitData(data1)
        split2 = splitData(data2)
        split3 = splitData(data2)
        
        learner1.train(split1['trainingdata'].rtable[:numEls, :], split1['traininglabels'][:numEls])
        learner2.train(split2['trainingdata'].rtable[:numEls, :], split2['traininglabels'][:numEls])
        
        stats1 = learner1.test(split2['testingdata'].rtable, split2['testinglabels'])
        tracker1.addStats(stats1['precision'], stats1['recall'], stats1['rc'])
        
        stats2 = learner2.test(split3['testingdata'].rtable, split3['testinglabels'])
        tracker2.addStats(stats2['precision'], stats2['recall'], stats2['rc'])
    
    losses = dict()
    
    losses['f1'] = tracker1.f1s/tracker2.f1s
    losses['area'] = tracker1.areas/tracker2.areas
    
        
    return losses
        
def output_crosstrain(specieslist, numEls = 1000, numtries = 1000):
    ctf1s = np.zeros((len(specieslist), len(specieslist)))
    ctf1s_pv = np.zeros((len(specieslist), len(specieslist)))
    
    ctars = np.zeros((len(specieslist), len(specieslist)))
    ctars_pv = np.zeros((len(specieslist), len(specieslist)))
    
    for s1 in range(len(specieslist)):
        for s2 in range(len(specieslist)): 
            losses = cross_train(specieslist[s1], specieslist[s2], numEls, numtries)
            ctf1s[s1,s2] = np.mean(losses['f1'])
            ctf1s_pv[s1, s2] = stats.ttest_1samp(losses['f1'], 1)[1]

            ctars[s1, s2] = np.mean(losses['area'])
            ctars_pv[s1, s2] = stats.ttest_1samp(losses['area'], 1)[1]
    
    np.savetxt('F1_CT.csv', ctf1s, delimiter='\t')
    np.savetxt('F1_CT_PV.csv', ctf1s_pv, delimiter='\t')
    
    np.savetxt('AR_CT.csv', ctars, delimiter='\t')
    np.savetxt('AR_CT_PV.csv', ctars_pv, delimiter='\t')