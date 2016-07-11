'''
myConstants.py --- various constants used by the program
'''

dbname = 'prova3'                                                                               # Name of the mySql database  
dbuser = 'gont'                                                                                 # Username to access the database
dbpasswd = 'gont'                                                                               # Password to access the database
ontodir = 'datasets/ontology/'                                                                  # Directory containing the ontology and slim ontology files
godir = ontodir+'go_daily-termdb-tables'                                                        # The subdirectory containing the ontology 
slimfile = ontodir+'goslim_generic.obo'                                                         # The slim file

allspecies = {'hsap', 'ecoli', 'mus', 'scere'}                                                  # Names of the species to study

annotsfiles = dict()                                                                            # Gene association files 
annotsfiles['hsap']='datasets/hsap/annotations/gene_association.goa_ref_human'
annotsfiles['ecoli']='datasets/ecoli/annotations/gene_association.ecocyc'
annotsfiles['mus']='datasets/mus/annotations/gene_association.mgi'
annotsfiles['scere']='datasets/scere/annotations/gene_association.sgd'


crossreffiles = dict()                                                                          # Cross-reference files (if proteins are not identified through uniprot IDs in the gene association files, otherwise NULL)
crossreffiles['hsap']=None
crossreffiles['ecoli']=None
crossreffiles['mus'] = 'datasets/mus/interactions/cr_mus.csv'
crossreffiles['scere'] = 'datasets/scere/interactions/cr_scere.csv'


coreintfiles = dict()                                                                           # Core (that is, trusted) interaction files 
coreintfiles['hsap'] = 'datasets/hsap/interactions/Hsapi20150701CR.txt'
coreintfiles['ecoli'] = 'datasets/ecoli/interactions/Ecoli20150701CR.txt'
coreintfiles['mus'] = 'datasets/mus/interactions/Mmusc20150701CR.txt'
coreintfiles['scere'] = 'datasets/scere/interactions/Scere20150701CR.txt'


fullintfiles = dict()                                                                           # Full (that is, also untrusted) interaction files
fullintfiles['hsap'] = 'datasets/hsap/interactions/Hsapi20150701.txt'
fullintfiles['ecoli'] = 'datasets/ecoli/interactions/Ecoli20150701.txt'
fullintfiles['mus'] = 'datasets/mus/interactions/Mmusc20150701.txt'
fullintfiles['scere'] = 'datasets/scere/interactions/Scere20150701.txt'


simfiles = dict()                                                                              # Where to store the interaction data for interacting pairs (core only)
simfiles['hsap'] = 'generated_files/pickled/hsap_idata.pickle'
simfiles['ecoli'] = 'generated_files/pickled/ecoli_idata.pickle'
simfiles['mus'] = 'generated_files/pickled/mus_idata.pickle'
simfiles['scere'] = 'generated_files/pickled/scere_idata.pickle'

simfilesfull = dict()                                                                          # Where to store the interaction data for interacting pairs (all full)
simfilesfull['hsap'] = 'generated_files/pickled/hsap_ifdata.pickle'
simfilesfull['ecoli'] = 'generated_files/pickled/ecoli_ifdata.pickle'
simfilesfull['mus'] = 'generated_files/pickled/mus_ifdata.pickle'
simfilesfull['scere'] = 'generated_files/pickled/scere_ifdata.pickle'


nonsimfiles = dict()                                                                           # Where to  store the interaction data for non-interacting pairs
nonsimfiles['hsap'] = 'generated_files/pickled/hsap_nidata.pickle'
nonsimfiles['ecoli'] = 'generated_files/pickled/ecoli_nidata.pickle'
nonsimfiles['mus'] = 'generated_files/pickled/mus_nidata.pickle'
nonsimfiles['scere'] = 'generated_files/pickled/scere_nidata.pickle'

statfiles = dict()                                                                      
statfiles['hsap'] = 'generated_files/hsap_stats.csv'
statfiles['ecoli'] = 'generated_files/ecoli_stats.csv'
statfiles['mus'] = 'generated_files/mus_stats.csv'
statfiles['scere'] = 'generated_files/scere_stats.csv'

pvalfiles = dict()
pvalfiles['hsap'] = 'generated_files/hsap_pvals.csv'
pvalfiles['ecoli'] = 'generated_files/ecoli_pvals.csv'
pvalfiles['mus'] = 'generated_files/mus_pvals.csv'
pvalfiles['scere'] = 'generated_files/scere_pvals.csv'

crosstrain_f1_file = 'generated_files/crosstrain_f1.csv'
crosstrain_f1_pv_file = 'generated_files/crosstrain_f1_pv.csv'
crosstrain_roc_file = 'generated_files/crosstrain_roc.csv'
crosstrain_roc_pv_file = 'generated_files/crosstrain_roc_pv.csv'



