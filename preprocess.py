'''
Preprocess.py --- functions to load the ontology and the species data, and generate the interacting pairs
'''

import subprocess
from myConstants import *
import MySQLdb
import csv
from collections import defaultdict




def loadgo():
    '''
        Loads the gene ontology
    '''
    subprocess.call('cat ' + godir + '/*.sql | mysql ' '--user=' + dbuser + ' --password='+dbpasswd +' ' + dbname, shell=True)
    subprocess.call('mysqlimport -L --user=' + dbuser + ' --password='+dbpasswd + ' ' + dbname + ' ' + godir + '/*.txt', shell = True)
    

def loadslims():
    '''
        Loads the slim sub-ontology
    '''
    slhandler = open(slimfile, 'r')
    db = MySQLdb.connect('127.0.0.1', dbuser, dbpasswd, dbname)
    cursor = db.cursor()
    cursor.execute('CREATE TABLE generic_slim (acc varchar(100), go_id int, slim_id int NOT NULL AUTO_INCREMENT, PRIMARY KEY(slim_id))')

    for line in slhandler: 
            line = line.strip()
            key, _, val = line.partition(":")     # Statements are of the form key:value  
            if key == 'id':                         # id of the current term
                val = val.strip()
                if val[:2]=='GO':
                    cursor.execute('INSERT INTO generic_slim VALUES (%s, NULL, NULL)', {val})
    
    cursor.execute('UPDATE generic_slim SET go_id = (SELECT id FROM term WHERE term.acc = generic_slim.acc)')
    db.commit()
    cursor.close()
    
                    
def loadannots(species):
    '''
        Loads the protein association data for the given species
    '''
    species_assoc = species + '_assoc'
    annfile = open(annotsfiles[species], 'r')
    db = MySQLdb.connect('127.0.0.1', dbuser, dbpasswd, dbname)
    cursor = db.cursor()
    cursor.execute('CREATE TABLE '+species_assoc +' (prod varchar(100), acc varchar(100), go_id varchar(100), PRIMARY KEY(prod, go_id))')

    for line in annfile: 
        if line[0]!='!': 
            columns = line.split('\t')
            qualifier = set(columns[3].strip().split('|'))
            if 'NOT' not in qualifier: 
                acc = columns[4].strip()
                prod_ID = columns[1]
                cursor.execute('INSERT IGNORE INTO '+species_assoc+' VALUES (%s, %s, (SELECT id FROM term WHERE term.acc=%s))', [prod_ID, acc, acc])
    db.commit()
    cursor.close()
                
def computeics(species):
    '''
        Computes information contents (for the given species)
    '''
    species_assoc = species+'_assoc'
    species_term_count = species + '_term_count'
    species_cumul_count = species + '_cumul_count'
    species_type_count = species + '_type_count'
    species_ic = species + '_ic'
    db = MySQLdb.connect('127.0.0.1', dbuser, dbpasswd, dbname)
    cursor = db.cursor()
    cursor.execute('create table '+species_term_count+' (term_id int not null, count int, primary key(term_id))')
    cursor.execute('insert into ' + species_term_count + ' (select term.id, count(*) from ' + species_assoc+ ' inner join term on go_id = term.id group by term.id)')
    cursor.execute('insert into ' + species_term_count + ' (select id, 0 from term where id not in (select distinct term_id from ' + species_term_count+ '))')
    cursor.execute('create table ' + species_cumul_count +' like ' +  species_term_count)
    cursor.execute('insert into ' + species_cumul_count + ' (select term1_id, sum(tc.count) from (select distinct term1_id, term2_id from graph_path\
                    inner join term as rel on relationship_type_id = rel.id) as gp inner join ' + species_term_count + ' as tc  on tc.term_id = term2_id group by term1_id);')
    cursor.execute('create table ' + species_type_count + ' (type_id int, name varchar(100), count int)')
    cursor.execute('insert into ' + species_type_count + ' (select tt.id, name, count from (select distinct term.id, term.name from term where term.name in\
                    (select distinct term_type from term)) as tt inner join ' + species_cumul_count + ' on tt.id = '+species_cumul_count+'.term_id);')
    
    cursor.execute('create table ' + species_ic + ' (term_id int not null, ic double, primary key(term_id));')
    cursor.execute('insert into ' + species_ic + ' (select term.id, -LOG2('+ species_cumul_count+'.count/' + species_type_count+'.count)/LOG2('+species_type_count+'.count) from term\
                    inner join ' + species_cumul_count+' on term.id = ' + species_cumul_count+'.term_id inner join '+ 
                    species_type_count+' on term.term_type = ' + species_type_count+'.name);')
    db.commit()
    cursor.close()
    
def loadcrfile(species):
    '''
        Loads the cross-reference file (for the given species)
    '''
    myfile = open(crossreffiles[species], 'r')
    myreader = csv.reader(myfile, delimiter='\t')
    cr = defaultdict(list)
    for l in myreader: 
        #print(str(l))
        if len(l)>1:
            cr[l[0]] += [l[1]]
    return cr
    
def loadintfile(species, intfile, tablename):
    '''
        Loads an interaction file (for a given species). intfile = interaction file name, tablename = name of corresponding mySql table
    '''
    if crossreffiles[species]!=None:
        crossreference = loadcrfile(species)
    else:
        crossreference = None
        
    species_assoc = species + '_assoc'
    myfile = open(intfile, 'r')
    db = MySQLdb.connect('127.0.0.1', dbuser, dbpasswd, dbname)
    cursor = db.cursor()
    cursor.execute('create table ' + tablename +  ' (prod1 varchar(100), prod2 varchar(100), primary key(prod1, prod2))')
    myreader = csv.reader(myfile, delimiter='\t')
    
    numnotfound = 0
    for l in myreader: 
        prod1 = l[0].split('|')[-1]
        prod2 = l[1].split('|')[-1]
        if prod1[:10]=='uniprotkb:' and prod2[:10]=='uniprotkb:':
            if crossreference == None: 
                cursor.execute('insert ignore into ' + tablename + ' values (%s, %s)', [prod1[10:], prod2[10:]])
            else: 
                for p1 in crossreference[prod1[10:]]:
                    if len(p1)==0:
                        #print('Uniprot not found: ' + prod1[10:])
                        numnotfound += 1
                    for p2 in crossreference[prod2[10:]]:
                        if len(p2)==0:
                            #print('Uniprot not found: ' + prod2[10:])
                            numnotfound += 1
                        cursor.execute('insert ignore into ' + tablename + ' values (%s, %s)', [p1, p2])
                
                

    numdel = cursor.execute('delete from ' + tablename + ' where prod1 not in (select distinct prod from ' + species_assoc + ') or prod2 not in (select prod from '+species_assoc+')')
    cursor.execute('select count(*) from ' + tablename)
    num_interacting = cursor.fetchall()[0][0]
    db.commit()
    cursor.close()
    #print('Uniprots not found in translation: ' + str(numnotfound) + ', not found in association file: ' + str(numdel))
    return num_interacting
    
def loadinteracting(species):
    '''
        Loads interaction data for the given species, generates non-interacting pairs
    '''
    species_int= species+'_int'
    species_intfull = species+'_int_full'
    species_notint = species+'_notint'
    species_intprods = species+'_intprods'
    
    loadintfile(species, coreintfiles[species], species_int)
    loadintfile(species, fullintfiles[species], species_intfull)
    
    db = MySQLdb.connect('127.0.0.1', dbuser, dbpasswd, dbname)
    cursor = db.cursor()
    cursor.execute('create table ' + species_intprods + ' (prod varchar(100), primary key(prod))')
    cursor.execute('insert ignore into ' + species_intprods + ' (select prod1 from ' + species_int+');')
    cursor.execute('insert ignore into ' + species_intprods + ' (select prod2 from ' + species_int+');')

    cursor.execute('create table ' + species_notint + ' like ' +species_int)
    for _ in range(3):
        cursor.execute('set @counter = 0;')
        cursor.execute('insert ignore into ' + species_notint + ' select prodA, prodB from (select prod as prodA, floor(rand()*(select count(*) from ' +
                       species_intprods + '))+1 as idxA from ' + species_intprods + ') as protA inner join (select prod as prodB, (@counter := @counter+1) as idxB from ' + 
                       species_intprods+ ') as protB on idxA = idxB;');
    cursor.execute('delete from ' + species_notint + ' where (prod1, prod2) in (select * from ' + species_intfull+')')
    cursor.execute('delete from ' + species_notint+' where (prod2, prod1) in (select * from '+species_intfull+ ')')
    db.commit()
    cursor.close()

    
def loadspecies(species):
    '''
        Loads all data pertaining to a given species
    '''
    print(species + ': loading annotations...')
    loadannots(species)
    print(species + ': computing information contents...')
    computeics(species)
    print(species + ': loading interaction data...')
    loadinteracting(species)
    
def loadall():
    print('Loading gene ontology...')
    loadgo()
    print('Loading slim sub-ontology...')
    loadslims()
    for species in allspecies: 
        loadspecies(species)
