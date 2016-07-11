INSTALLATION NOTES

Because of upload size limitations, the data has been split into parts and separated from the code and the generated data has not been included.

To run the code, first set the name of the mySql database (which must already exist) and the user name and password for it in the file myConstants.py. Then build the database and perform the preprocessing calling the function prepareAll() in the file logsim.py.

To test the performance of the similarity over a species, just execute testSpecies from logsim with respect to the corresponding species abbreviation (species available at the moment: 'scere' = S. cerevisiae, 'hsap' = H. sapiens, 'ecoli' = E. coli, 'mus' = M. musculus). The result will be put in the generated_files folder, in the two files <species>_stats.csv and <species>_pvals.csv. The former contains the means and the standard deviations of a few performance measures (precision, recall, F1, area under ROC) for the new measure and for the Resnik, Lin, Jiang-Conrath and simGIC measures (in this order)
evaluated first all wrt the biological process sub-ontology, then all wrt the molecular function sub-ontology, and finally all wrt the cellular component sub-ontology. The second file contains the t-statistics and the p-values obtained by performing (for every choice of evaluation metric) dependent t-tests with paired samples between the new measure and the best of the non-adaptive measures. 

The function returns a simTester objects - see the file simTester.py for its complete definition. Two of its methods which hold particular interest are: 

a. plotrcs(root): plots the ROC curves relative to the adaptive similarity and to the non-adaptive ones, these latter calculated over the sub-ontology root = 'BIO', 'MOL' or 'CEL'.

b. plot_weights(): plots the mean values of the weights that the adaptive similarity learned for the GO-slim terms. 

NOTE: for the default number of tries (1000) per test, this operation may last quite a while. The program will print an output at every iteration until 100 iterations are reached; then it will print an output every 100 iterations. 

3) To test the performance of cross-training, run crossTrain (also declared in logSim). NOTE: it will take a LONG time to finish, for the default value of numEl (number of elements used for the training sets) of 1000. The output will be stored inside the generated_files folder, in the files crosstrain_<f1/roc>.pv for the gain/losses and in the files crosstrain_<f1/roc>_pv.csv for the corresponding p-values against the hypothesis gain=1 (that is, no significant difference between training and cross-training).

4) It is also possible to train directly the adaptive similarity measure, which is declared in simLearner. 

Example: 
> hdata = simLearner.getData('hsap')
> hsplit = simLearner.splitData(hdata)
> sl = simLearner.simLearner()
> sl.train(hsplit['trainingdata'].rtable, hsplit['traininglabels'])
> sl.test(hsplit['testingdata'].rtable, hsplit['testinglabels'])

Happy testing!
