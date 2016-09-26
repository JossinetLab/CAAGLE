#!/usr/bin/env python

import sys, os, urllib, urllib2
from zipfile import ZipFile
from pymongo import MongoClient
from bs4 import BeautifulSoup


def download_data():
	dirname = os.path.abspath(os.path.dirname(__file__))+'/../data'
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	### BIOGRID DATA ###
	url = "http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.138/BIOGRID-ORGANISM-3.4.138.tab2.zip"
	biogrid_path = dirname+'/BioGrid'
	if not os.path.exists(biogrid_path):
		os.makedirs(biogrid_path)
	zip_response = urllib.urlopen(url)
	tempory_file = url.split('/')[-1]
	with open("%s/%s"%(biogrid_path, tempory_file), "wb") as fh:
		fh.write(zip_response.read())
	zip_response.close()
	zf = ZipFile("%s/%s"%(biogrid_path, tempory_file))
	zf.extractall(path=biogrid_path)
	zf.close()
	os.remove("%s/%s"%(biogrid_path, tempory_file))
	biogrid_file_name = None
	for file_name in os.listdir(biogrid_path):#we upload a directory
		if not 'Saccharomyces_cerevisiae_S288c' in file_name:
			os.remove("%s/%s"%(biogrid_path, file_name))
		else:
			biogrid_file_name = file_name

	### SGD DATA ###
	url = "http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
	sgd_path = dirname+'/SGD'
	if not os.path.exists(sgd_path):
		os.makedirs(sgd_path)
	response = urllib.urlopen(url)
	outfile_content = str(response.read())
	response.close()
	sgd_file_name = url.split('/')[-1]
	with open("%s/%s"%(sgd_path, sgd_file_name), 'w') as fh:
		fh.write(outfile_content)

	return biogrid_path+"/"+biogrid_file_name, sgd_path+"/"+sgd_file_name

def parse_biogrid_data(data_path):
	genetic_interactions = {}
	physical_interactions = {}
	with open(data_path, 'r') as fh:
		lines = fh.readlines()
		for line in lines:
			if not line.startswith('#'):
				tokens = line.strip().split('\t')
				interactorA = tokens[5]#BioGrid and SGD Systematic name
				interactorB = tokens[6]#BioGrid and SGD Systematic name
				interaction_type = tokens[12]#'genetic' or 'physical'
				interactorA_organism = tokens[15]
				interactorB_organism = tokens[16]
				if interactorA_organism == '559292' and interactorB_organism == '559292':#NCBI taxonomy ID of 'Saccharomyces_cerevisiae_S288c' 
					biogrid_dict = genetic_interactions if interaction_type == 'genetic' else physical_interactions
					intB_list = biogrid_dict.get(interactorA)
					if not intB_list:
						biogrid_dict[interactorA] =	[interactorB]
					else:
						if interactorB not in intB_list:#we remove redundancy
							intB_list.append(interactorB)

					intA_list = biogrid_dict.get(interactorB)
					if not intA_list:
						biogrid_dict[interactorB] =	[interactorA]
					else:
						if interactorA not in intA_list:#we remove redundancy
							intA_list.append(interactorA)

	return genetic_interactions, physical_interactions

def parse_sgd_data(data_path):
	correlation_id_name = {}
	with open(data_path, 'r') as fh:
		lines = fh.readlines()
		for line in lines:
			tokens = line.strip().split('\t')
			correlation_id_name[tokens[0].strip()] = tokens[3].strip()
	return correlation_id_name	

def search_interactants(genetic_interactions, physical_interactions, correlation_id_name):
	"""
    For each C. glabrata/Cagl annotation:
	- retrieve the SGD link of the S. cerevisiae/Sace ortholog (if it exists ; several? currently NOT)
	- convert the SGD ID, contained in the SGD link, into SGD systematic name
	- get Sace GENETIC and PHYSICAL BioGrid interactors 
		For each Sace interactor:
		- convert the SGD systematic name of this interactor into SGD ID
		- from the SGD ID, contained in the SGD link, search the entire Cagl MongoDB for the absence or the presence of one or several Cagl orthologs
		
		If no Cagl ortholog corresponds to the Sace interactor => keep the SGD link of the Sace interactor
		If one or several Cagl ortholog(s) => keep the MongoDB ID of the Cagl annotation(s) ; for example: "5774cd8f9ae2c726d247b7cf@Candida_glabrata_CBS_138"


    If interactors, update the C. glabrata MongoDB table named 'annotations':
    	'interactors':
    			'source': list of source(s) "db:biogrid:sgd_id"
    			'genetic': list of the Cagl ortholog ID of the Sace GENETIC interactors or, if no Cagl ortholog exists, the SGD link of the Sace interactor
    			'physical': list of the Cagl ortholog ID of the Sace PHYSICAL interactors or, if no Cagl ortholog exists, the SGD link of the Sace interactor

    """
	client = MongoClient()
	db = client['Candida_glabrata_CBS_138']
	
	sgd_url ="http://www.yeastgenome.org/cgi-bin/locus.pl?dbid="

	ids_to_update = []
	for annotation in db['annotations'].find({'orthologs_in_non_CGD_species':{'$regex':"S. cerevisiae:"}}, no_cursor_timeout = True):
		if annotation.get('interactors'):
			print "%s, interactors already stored"%annotation['locus_tag']
		else:
			print annotation['locus_tag']
			genetic_cagl_interactors = []
			physical_cagl_interactors = []
			source = []
			for ortholog in annotation['orthologs_in_non_CGD_species']:
				if ortholog.split(':')[0] == 'S. cerevisiae':
					sgd_id = ':'.join(ortholog.split(':')[1:]).split('dbid=')[1]
					sgd_systematic_name = correlation_id_name.get(sgd_id)
					if sgd_systematic_name:
						
						### Search C. glabrata orthologs of S. cerevisiae GENETIC interactors ###
						genetic_interactors = genetic_interactions.get(sgd_systematic_name)
						if genetic_interactors:
							for genetic_interactor in genetic_interactors:#genetic_interactor = SGD Standard Name
								for genetic_item in correlation_id_name.iteritems():
									if genetic_item[1] == genetic_interactor:
										if db['annotations'].find({'orthologs_in_non_CGD_species': {'$regex': "dbid=%s"%genetic_item[0]}}).count():
											for genetic_annotation in db['annotations'].find({'orthologs_in_non_CGD_species': {'$regex': "dbid=%s"%genetic_item[0]}}, no_cursor_timeout = True):
												genetic_cagl_interactors.append("%s@Candida_glabrata_CBS_138"%genetic_annotation['_id'])
										else:# NO Cagl ortholog of the Sace genetic interactor 
											genetic_cagl_interactors.append(sgd_url+genetic_item[0])
						
						### Search C. glabrata orthologs of S. cerevisiae PHYSICAL interactors ###
						physical_interactors = physical_interactions.get(sgd_systematic_name)
						if physical_interactors:
							for physical_interactor in physical_interactors:#physical_interactor = SGD Standard Name
								for physical_item in correlation_id_name.iteritems():
									if physical_item[1] == physical_interactor:
										if db['annotations'].find({'orthologs_in_non_CGD_species': {'$regex': "dbid=%s"%physical_item[0]}}).count():
											for physical_annotation in db['annotations'].find({'orthologs_in_non_CGD_species': {'$regex': "dbid=%s"%physical_item[0]}}, no_cursor_timeout = True):
												physical_cagl_interactors.append("%s@Candida_glabrata_CBS_138"%physical_annotation['_id'])
										else:# NO Cagl ortholog of the Sace physical interactor
											physical_cagl_interactors.append(sgd_url+physical_item[0])
						
						if genetic_interactors or physical_interactors: 
							source.append("db:biogrid:%s"%sgd_id)

			if genetic_cagl_interactors or physical_cagl_interactors:
				interactors = {'source': source}
				if genetic_cagl_interactors: 
					interactors ['genetic'] = list(set(genetic_cagl_interactors))#we remove redundancy
				if physical_cagl_interactors:
					interactors ['physical'] = list(set(physical_cagl_interactors))#we remove redundancy
				ids_to_update.append((annotation['_id'], interactors))
	
	### Update C. glabrata MongoDB ###
	for _id in ids_to_update: 
		db['annotations'].find_one_and_update({'_id':_id[0]},{'$set':{'interactors':_id[1]}}, upsert=False)

	client.close()

if __name__ == '__main__':
	biogrid_path, sgd_path = download_data()
	genetic_interactions, physical_interactions = parse_biogrid_data(biogrid_path)
	#genetic_interactions, physical_interactions = parse_biogrid_data("/Users/laurence/CAAGLE/scripts/../data/BioGrid/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.138.tab2.txt")
	correlation_id_name = parse_sgd_data(sgd_path)
	#correlation_id_name = parse_sgd_data("/Users/laurence/CAAGLE/scripts/../data/SGD/SGD_features.tab")
	
	# search_interactants(parse_biogrid_data(biogrid_path), parse_sgd_data(sgd_path)) #FUNCTION search_interactants() REQUEST 3 ARGUMENTS !!!!!!!!! genetic_interactions, physical_interactions, correlation_id_name
	search_interactants(genetic_interactions, physical_interactions, correlation_id_name)
	
