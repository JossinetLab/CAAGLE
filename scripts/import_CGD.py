#!/usr/bin/env python

"""
COMMAND-LINE
------------
You can run the program with this command-line:
./import_CGD.py

For example:
./import_CGD.py -locus 1,10

ARGUMENTS
---------
-locus: list of Integers. Each Integer corresponds to the position of the locus tag in the file named locus_tag.txt (by default: None ; all locus tags will be taken into account).
    For example, with the argument -locus equal to 1,10 the function import_cgd() will only process the first ten locus tags of the file locus_tag.txt. 
"""

import os, sys, re, urllib, urllib2, gzip
from pyrna.parsers import parse_fasta, to_fasta
from pyrna.features import Molecule, DNA
from pyrna.computations import Tool
from pymongo import MongoClient
from bson import ObjectId
from bs4 import BeautifulSoup
import socket
from pandas import DataFrame

def dump_data():
    """    
    This function recovers all the data necessary to construct the CAGGLE DB:
        - from the databases GRYC (for the Nakaseomyces, http://www.igenolevures.org) and
        - from the CGD (for all the other species, http://www.candidagenome.org/) 
    It will produce one directory per species. In each directory, it will dump the genomic sequences in a file named sequences.fasta. 
    For the Nakaseomyces, it will also produce a file named annotations.embl.
    And for the species C. glabrata, this script will also produce a file named locus_tag.txt. This file will be used for the next step (see function import_cgd())
    """
    dirname = os.path.abspath(os.path.dirname(__file__))+'/../data'
    if not os.path.exists(dirname):
        # we create the data directory
        os.makedirs(dirname)
    # the species we're interested in
    species_dict = {
                    'Candida_glabrata_CBS_138': ['C_glabrata_CBS138'],
                    'Nakaseomyces_bacillisporus_CBS_7720': ['NABA0S', '01', 37],
                    'Nakaseomyces_bracarensis_CBS_10154': ['CABR0S', '01', 40],
                    'Nakaseomyces_castellii_CBS_4332': ['CACA0S', '01', 36],
                    'Nakaseomyces_delphensis_CBS_2170': ['NADE0S', '01', 30],
                    'Nakaseomyces_nivariensis_CBS_9983': ['CANI0S', '01', 29],
                    'Candida_albicans_SC5314': ['C_albicans_SC5314'],
                    'Candida_parapsilosis_CDC317': ['C_parapsilosis_CDC317'],
                    'Candida_dubliniensis_CD36': ['C_dubliniensis_CD36']
                    }
    candida = ['Candida_glabrata_CBS_138', 'Candida_albicans_SC5314', 'Candida_parapsilosis_CDC317', 'Candida_dubliniensis_CD36']
                  
    for item in species_dict.iteritems():
        species = item[0]
        data_path = dirname+'/'+species
        file_ids = item[1]
        if species not in candida:# we generate the list of chromosome sequence IDs for the Nakaseomyces species
            sequence_ids = []
            for n in range(file_ids[2]):
                sequence_ids.append(file_ids[0]+"%02u"%(int(file_ids[1])+n))
            file_ids = sequence_ids

        if not os.path.exists(data_path):
            # we create the species directory
            os.makedirs(data_path)
        # we recover and dump the data in this directory
        if species in candida:
            if species == 'Candida_albicans_SC5314':
                url = "http://www.candidagenome.org/download/sequence/%s/Assembly22/current/%s_A22_current_chromosomes.fasta.gz"%(file_ids[0], file_ids[0])
            else:
                url = "http://www.candidagenome.org/download/sequence/%s/current/%s_current_chromosomes.fasta.gz"%(file_ids[0], file_ids[0])
            response = urllib.urlopen(url)
            outfile_content = str(response.read())
            response.close()
            gz_file_name = url.split('/')[-1]
            with open("%s/%s"%(data_path, gz_file_name), 'w') as fh:
                fh.write(outfile_content)
            with gzip.open("%s/%s"%(data_path, gz_file_name), 'rb') as fh:
                file_content = fh.read()
            os.remove("%s/%s"%(data_path, gz_file_name))
            with open("%s/%s"%(data_path, 'sequences.fasta'), 'w') as fh:
                fh.write(file_content)

            if species == 'Candida_glabrata_CBS_138':
                outfile_content = ""
                for string in ['orf_genomic', 'other_features_genomic']:
                    url = "http://www.candidagenome.org/download/sequence/%s/current/%s_current_%s.fasta.gz"%(file_ids[0], file_ids[0], string)
                    response = urllib.urlopen(url)
                    outfile_content += str(response.read())
                    response.close()
                    gz_file_name = url.split('/')[-1]
                    with open("%s/%s"%(data_path, gz_file_name), 'w') as fh:
                        fh.write(outfile_content)
                    with gzip.open("%s/%s"%(data_path, gz_file_name), 'rb') as fh:
                        file_content = fh.read()
                    os.remove("%s/%s"%(data_path, gz_file_name))
                    with open("%s/%s"%(data_path, 'locus_tags.txt'), 'w') as fh:
                        dna_molecules = parse_fasta(file_content, 'DNA')
                        for d in dna_molecules:
                            fh.write("%s\n"%d.name.split()[0])            
        else:# Nakaseomyces
            match = re.compile('^([A-Za-z]+_[A-Za-z]+)_(.+)$').search(species)
            if match:
                for file_type in [['Fasta', 'sequences.fasta', 'fsa'], ['EMBL', 'annotations.embl', 'embl']]:
                    outfile_content = ""
                    for file_id in file_ids:
                        url = "http://gryc.inra.fr/storedData/download/%s/%s/%s/%s.%s"%(match.group(1), match.group(2), file_type[0], file_id, file_type[2])
                        response = urllib.urlopen(url)
                        outfile_content += str(response.read())
                        response.close()
                        with open("%s/%s"%(data_path, file_type[1]), 'w') as fh:
                            fh.write(outfile_content)

    print "The file downloading is done."

def parse_embl(embl_data, feature_keys=None):
    """
    This function retrieves all features and all its associated qualifiers in a text at the EMBL format. It returns a dataframe.

    Arguments:
    ---------
    -feature_keys: list of feature keys (default: None ; no sorting by feature key).
    Only the features having the keys given in argument will be returned by the function parse_embl().
    for example:
        feature_keys=['CDS','tRNA', 'ncRNA', 'rRNA']
    """
    all_cds = []
    lines = embl_data.strip().split('\n')
    if not lines[-1].strip() == '//':
        raise Exception("Uncomplete file")
    
    genome = ''
    location = ''
    feature = {}
    features = []
    qualifiers = []
    pattern1 = re.compile('^FT\s+([a-zA-Z_]+)\s+[a-z(]*([0-9.,><]+)[\s)]*$')
    pattern2 = re.compile('^FT\s+/([a-zA-Z_]+)=?(.+)$')#we include the 2nd double quote if it exists ; reg exp NOT correct: '^FT\s+/([a-zA-Z_]+)=(.+)$'
    pattern3 = re.compile('^FT\s+(.+)$')#a space is not at the beginning of the text but at the end (text is left aligned) ; we include the 2nd double quote if it exists
    for line in lines:
        if line.startswith('ID'):#we get all ID lines with the scaffold or chromosome name
            genome = line.split()[1].replace(';', '')
        ###  WE GET ALL FEATURES AND THEIR ASSOCIATED LOCATION AND QUALIFIERS IN THE FT LINES ###
        elif line.startswith('FT'):
            match1 = pattern1.search(line)
            if match1:#we get all FT lines with the feature type and feature genomic location
                if feature:
                    location = location.replace('>','')
                    location = location.replace('<','')
                    location = location.replace(')','')
                    positions = []
                    for exon in location.split(','):
                        if '..' in exon:
                            positions.append([int(exon.split('..')[0]), int(exon.split('..')[1])])
                        else:#feature with a unique coordinate! for example: locus_tag="BN123_NABA0s31e00254r"
                            positions.append([int(exon), int(exon)])
                    feature['location'] = positions#list of lists
                    if qualifiers:
                        for qualifier in qualifiers:
                            feature[qualifier[0]] = qualifier[1].replace('"', '').strip() if qualifier[0] != 'estimated_length' else int(qualifier[1])
                        qualifiers = []
                    features.append(feature)
                    feature = {}
                feature['feature_type'] = match1.group(1)
                location = match1.group(2)
                feature['strand'] = '-' if 'complement' in line else '+'
                feature['genome'] = genome.strip()
            else:
                match2 = pattern2.search(line)
                if match2:#we get all FT lines with a feature qualifier (for example: /locus_tag=)
                    qualifiers.append([match2.group(1), match2.group(2)])#we append qualifier and description
                else:#we get the other FT lines (without feature type nor feature qualifier)
                    match3 = pattern3.search(line)
                    if match3:
                        if qualifiers:#we complete the description of a qualifier 
                            current_qualifier = qualifiers.pop()
                            description = current_qualifier[1]+match3.group(1) if current_qualifier[0] == 'translation' else current_qualifier[1]+' '+match3.group(1) 
                            qualifiers.append([current_qualifier[0], description])
                        else:#we complete the location of a feature ; for example: locus_tag="BN121_CACA0s22e00110g"
                            location += match3.group(1)
        elif line.startswith('//') and feature:#we get the last feature of a scaffold or chromosome
            location = location.replace('>','')
            location = location.replace('<','')
            location = location.replace(')','')
            positions = []
            for exon in location.split(','):
                if '..' in exon:
                    positions.append([int(exon.split('..')[0]), int(exon.split('..')[1])])
                else:
                    positions.append([int(exon), int(exon)])
            feature['location'] = positions
            if qualifiers:
                for qualifier in qualifiers:
                    feature[qualifier[0]] = qualifier[1].replace('"', '').strip() if qualifier[0] != 'estimated_length' else int(qualifier[1])
                qualifiers = []
            features.append(feature)
            feature = {}
    
    ### FEATURE KEY SORTING ###
    if feature_keys:
        sorted_features = []
        for f in features:
            for f_type in feature_keys:
                if f['feature_type'] == f_type:
                    sorted_features.append(f)
        return DataFrame(sorted_features)#sorted_features is a list of dict
    else:
        return DataFrame(features)#features is a list of dict

def create_databases():
    """
    CGD: Candida Genome Database http://www.candidagenome.org/
    GRYC: Genome Resources for Yeast Chromosomes http://gryc.inra.fr/

    Step 1. For the four CGD Candida species (C. glabrata CBS138, C. albicans SC5314, C. dubliniensis CD36, C. parapsilosis CDC317) and the five GRYC Nakaseomyces species (C. bracarensis, C. castellii, C. nivariensis, N. bacillisporus, N. delphensis),
    creation of a MongoDB with the name of each species that contains a table named 'genomes' with the following fields:
        '_id'
        'name'
        'sequence'

    Step 2. For the five GRYC Nakaseomyces species,
    For the MongoDB of these species, creation of a table named 'annotations' with the following fields:
        '_id'
        'source': "db:gryc:species_name"
        'feature_type': feature key ('CDS' or 'tRNA' or 'ncRNA' or 'rRNA') 
        'locus_tag'
        'description': gene annotation of GRYC
        'genome': "genome_id@genomes"
        'genomeName'
        'genomicStrand'
        'genomicPositions': {'CDS' or 'Noncoding_exon': [[exon_start, exon_end]]}
        'sequence': genomic sequence. If intron(s), genomic spliced sequence
        'translation': protein sequence if feature key is 'CDS'
        'ncRNA_class': if feature key is 'ncRNA'

    Comments:
    --------
        list of all existing feature keys: ['LTR', 'ncRNA', 'repeat_region', 'tRNA', 'gap', 'source', 'misc_feature', 'rRNA', 'CDS', 'gene']
        list of CDS qualifiers used: ['locus_tag','translation','genome','strand','location','note','product']
        list of RNA qualifiers used: ['locus_tag','genome','strand','location','gene','product','note']
    """
    client = MongoClient()

    dirname = os.path.abspath(os.path.dirname(__file__))+'/../data'

    pattern = re.compile('^BN\d{3,3}_(\w{4,4}\ds\d\d)e.+$')

    for species in os.listdir(dirname):
        if os.path.isdir("%s/%s"%(dirname, species)):
            if "sequences.fasta" in os.listdir(dirname+'/'+species):
                print "Creation of the table 'genomes' for the database: %s"%species
                with open("%s/%s/sequences.fasta"%(dirname, species)) as h:
                    content = h.read()
                    dna_molecules = parse_fasta(content, 'DNA')
                    for d in dna_molecules:
                        client[species]['genomes'].insert(
                            {
                                '_id': str(ObjectId()),
                                'name': d.name.split()[0],
                                'sequence': d.sequence
                            }
                        )
            if "annotations.embl" in os.listdir(dirname+'/'+species):
                print "Creation of the table 'annotations' for the database: %s"%species
                with open("%s/%s/annotations.embl"%(dirname, species)) as h:
                    content = h.read()
                    df = parse_embl(content, feature_keys=['CDS','tRNA', 'ncRNA', 'rRNA'])
                    if not df.empty:
                        for row in df.iterrows():
                            data_embl = row[1]    
                            infos = {
                                    '_id': str(ObjectId()),
                                    'source': "db:gryc:%s"%species,
                                    'feature_type': data_embl['feature_type'],#feature key
                                    'locus_tag': data_embl['locus_tag'],
                                    'genomeName': data_embl['genome'],#for examples: "NABA0S01", "CANI0S01"
                                    'genomicStrand': data_embl['strand'],#'+' or '-'
                                    'genomicPositions': {data_embl['feature_type']: data_embl['location']} if data_embl['feature_type'] == 'CDS' else {'Noncoding_exon': data_embl['location']}
                                    }
                            description = ''
                            if infos['feature_type'] == 'CDS':
                                infos['translation'] = data_embl['translation']#protein sequence always present
                                if data_embl.notnull().get_value('note'):
                                    description = data_embl['note']
                                else:
                                    if data_embl.notnull().get_value('product'):
                                        description = data_embl['product']
                            elif infos['feature_type'] == 'ncRNA':
                                if data_embl.notnull().get_value('product'):
                                    description = data_embl['product']
                                if data_embl.notnull().get_value('ncRNA_class'):
                                    infos['ncRNA_class'] = data_embl['ncRNA_class']
                            elif infos['feature_type'] == 'tRNA' or infos['feature_type'] == 'rRNA':
                                if data_embl.notnull().get_value('gene'):
                                    description = data_embl['gene']
                                else:
                                    if data_embl.notnull().get_value('product'):
                                        description = data_embl['product']     
                                    else:
                                        if data_embl.notnull().get_value('note'):
                                            description = data_embl['note']
                                            if ';' in description and infos['feature_type'] == 'tRNA':
                                                description = description.split(';')[-1].strip()
                            if description:
                                infos['description'] = description
                            genome = client[species]['genomes'].find_one({'name': infos['genomeName']})
                            dna = DNA(name = genome['name'], sequence = genome['sequence'])
                            if genome:
                                infos['genome'] = genome['_id']+'@genomes'
                                if len(data_embl['location']) == 1:#gene without intron
                                    infos['sequence'] = dna[data_embl['location'][0][0]-1:data_embl['location'][0][1]] if data_embl['strand'] == '+' else dna.get_complement()[data_embl['location'][0][0]-1:data_embl['location'][0][1]][::-1]
                                else:#gene with intron(s)
                                    spliced_sequence = ''
                                    for exon in data_embl['location']:
                                        spliced_sequence += dna[exon[0]-1:exon[1]]
                                    infos['sequence'] = spliced_sequence if data_embl['strand'] == '+' else DNA(name=data_embl['locus_tag'], sequence=spliced_sequence).get_complement()[::-1]
                            client[species]['annotations'].insert(infos)    
    client.close()
    print "The database creation is done."

def import_cgd(fh, locus_range):
    """
    From each gene of Candida glabrata CBS138, some informations are extracted from its HTML page on the site CGD (Candida Genome Database).

    Creation of a table named 'annotations' for the four CGD Candida species with the following fields:
        '_id'
        'locus_tag'
        'source': "db:cgd:locus_tag"
        'standard_name'
        'alias'
        'genomicStrand'
        'genome': "genome_id@genomes"
        'genomeName'
        'genomicPositions': dictionary of subfeatures (5' UTR, CDS, Intron, 3' UTR, Noncoding_exon, Centromere, Long_terminal_repeat or Repeat_region) and their chromosomal coordinates
        'sequence': genomic sequence. If intron(s), genomic spliced sequence. If CDS and UTR, genomic translated sequence.
        'feature_type'
        'description'
        'orthologs_in_candida_species': list of "annotation_id@species_name"
        'orthologs_in_non_CGD_species': list of HTTP links
        'go_annotations': dictionary of GO terms and their associated description
        'external_links'
        'translation': protein sequence if the annotation is a CDS
        '3Dstructure': list of PDB links
        'sace_ortholog': list with the systematic name and protein sequence of the S. cerevisiae ortholog ; useful subsequently to do a multiple sequence alignment (see the script search_orthologs.py)
    """
    client = MongoClient()

    db = client['Candida_glabrata_CBS_138']

    dirname = os.path.abspath(os.path.dirname(__file__))+'/../data'

    with open("%s/Candida_glabrata_CBS_138/locus_tags.txt"%dirname) as h:
        if locus_range:
            for locus_tag in h.readlines()[int(locus_range[0])-1:int(locus_range[1])]:
                parse_cgd_entry(client, db, locus_tag.strip(), fh)# parse_cgd_entry() is a recursive function
        else:
            for locus_tag in h.readlines():
                parse_cgd_entry(client, db, locus_tag.strip(), fh)

    client.close()

def parse_cgd_entry(client, db, locus_tag, fh):

    cgd_infos = db['annotations'].find_one({'locus_tag': locus_tag})

    if cgd_infos:
        print locus_tag, "already stored"
        return cgd_infos['_id']

    mongodatabases = {
                    'C. dubliniensis CD36': 'Candida_dubliniensis_CD36',
                    'C. parapsilosis CDC317': 'Candida_parapsilosis_CDC317',
                    'C. albicans SC5314': 'Candida_albicans_SC5314',
                    'C. glabrata CBS138': 'Candida_glabrata_CBS_138'
                    }
    go_types = {
                "Molecular Function": 'F',
                "Biological Process": 'P',
                "Cellular Component": 'C'
                }

    cgd_infos = {
                'source': 'db:cgd:%s'%locus_tag,
                '_id': str(ObjectId()),
                'locus_tag': locus_tag
                }

    if locus_tag.startswith('orf') and db.name == 'Candida_albicans_SC5314': # see the comment below at level of the condition: if text == "Systematic Name, Reference Strain"
        pass
    else:
        db['annotations'].insert(cgd_infos)

    url = "http://www.candidagenome.org/cgi-bin/locus.pl?locus=%s"%locus_tag.strip()
    print "Entry parsing: ", locus_tag

    candida_orthologs = []
    non_cgd_orthologs = []
    go_terms = {}
    external_links = []
    try:
        soup = BeautifulSoup(str(urllib2.urlopen(url, timeout = 60.0).read()))
        for a in soup.find_all(id="tabtext", href=True):
            if a.next_sibling == ' Protein':
                protein_soup = BeautifulSoup(str(urllib.urlopen(a.get('href')).read()))
                for td in protein_soup.find_all("td", {"bgcolor":"#FFCC33"}):
                    if td.get_text().strip() == 'in FASTA format)':
                        params = urllib.urlencode({'downloadFasta': '1', 'feature_name': td.find("form").find_next('input').find_next('input').get('value')})
                        response = urllib2.urlopen("http://www.candidagenome.org%s"%td.find("form").get('action'), params)
                        fasta_protein = response.read()
                        cgd_infos['translation'] = parse_fasta(fasta_protein, 'DNA')[0].sequence.replace('*', '')
        for td in soup.find_all("td", {"bgcolor":"#FFCC33", 'valign':'top'}):
            text = td.get_text().strip()
            if text == "Standard Name":
                cgd_infos['standard_name'] = td.next_sibling.next_sibling.get_text("|").split("|")[0]
            elif text == "Systematic Name, Reference Strain":
                """
                To eliminate redundancy with the locus tags in the Candida albicans MongoDB
                For example: Systematic Name 'C1_01530C_A' = Assembly 19/21 file_id 'orf19.3341'
                """
                if locus_tag.startswith('orf') and db.name == 'Candida_albicans_SC5314': #for example: 'orf19.2769'
                    locus_tag = td.next_sibling.next_sibling.get_text().split()[0]
                    cgd_infos['locus_tag'] = locus_tag #'orf19.2769' becomes 'C4_02340W_A'
                    cgd_infos['source'] = 'db:cgd:%s'%locus_tag #'db:cgd:orf19.2769' becomes 'db:cgd:C4_02340W_A' and cgd_infos['_id'], cgd_infos['translation'] doesn't change  
                    cgd_annotation = db['annotations'].find_one({'locus_tag': locus_tag})
                    if cgd_annotation:
                        print locus_tag, "already stored"
                        return cgd_annotation['_id']
                    else:#cgd_infos is retained to continue
                        db['annotations'].insert(cgd_infos)
                else:
                    pass
            elif text == "Feature Type":
                cgd_infos['feature_type'] = td.next_sibling.next_sibling.get_text().strip()
            elif text == "Description":
                description = td.next_sibling.next_sibling.get_text().split('Literature')[0].strip()
                if re.search(' \([0-9, ]*\)$', description):
                    cgd_infos['description'] = re.split(' \([0-9, ]*\)$', description)[0]
            elif text == "Alias":
                alias = []
                for string in re.split('\|[0-9|]*', td.next_sibling.next_sibling.get_text("|")):
                    alias += string.split(', ')
                while '' in alias:
                    alias.remove('')
                cgd_infos['alias'] = alias
            elif text == "Best hit(s) in non-CGD species" or text == "Ortholog(s) in non-CGD species":
                sace_protein = None
                pdb_links = []
                for i in td.next_sibling.next_sibling.find_all('i'):
                    species = i.get_text().strip()
                    http_link = i.find_next('a', href=True).get('href')
                    non_cgd_orthologs.append("%s:%s"%(species, http_link))
                    ### we get the PDB link(s) of the Sace ortholog 3D structure(s), if it/they exist(s) ###
                    if species == 'S. cerevisiae':
                        sgd_soup = BeautifulSoup(str(urllib.urlopen(http_link).read()))
                        if sgd_soup.find(id="overview").find('dl', class_="key-value").find('dt', text="Systematic Name"):
                            sgd_systematic_name = sgd_soup.find(id="overview").find('dl', class_="key-value").find('dt', text="Systematic Name").find_next('dd').get_text().strip()
                        if sgd_soup.find(id="resources"):
                            for a in sgd_soup.find(id="resources").find_all('a', href=True):
                                if a.get_text().strip() == 'UniProtKB':
                                    uniprot_soup = BeautifulSoup(str(urllib.urlopen(a.get('href')).read()))
                                    if uniprot_soup.find(id="structure"):
                                        if uniprot_soup.find(id="structure").find(id="section_x-ref_structure"):
                                            for a in uniprot_soup.find(id="structure").find(id="section_x-ref_structure").find_next('table', class_="databaseTable STRUCTURE").find_all('a', href=True):
                                                if a.get('class'):#a.get('class') returns a list or None
                                                    for c in a.get('class'):
                                                        if c == 'pdb':
                                                            pdb_links.append(a.get('href'))
                                    if uniprot_soup.find(id="sequences"):
                                        if uniprot_soup.find(id="sequences").find(class_="sequence-isoform-leftcol"):
                                            for a in uniprot_soup.find(id="sequences").find(class_="sequence-isoform-leftcol").find_all('a', href=True):
                                                if a.get('href').endswith(".fasta"):
                                                    sgd_protein_molecules = parse_fasta(str(urllib.urlopen("http://www.uniprot.org%s"%a.get('href')).read()), 'DNA')
                                                    sace_protein = (sgd_systematic_name, sgd_protein_molecules[0].sequence.replace('*', ''))#tuple(systematic name, protein sequence)
                                elif a.get_text().strip() == 'UniParc':
                                    uniparc_soup = BeautifulSoup(str(urllib.urlopen(a.get('href')).read()))
                                    if uniparc_soup.find(id="entrySequence"):
                                        fasta_tokens = uniparc_soup.find(id="entrySequence").get_text().strip().split('\n')
                                        if fasta_tokens[0].startswith('>'):
                                            sace_protein = (sgd_systematic_name, ''.join(fasta_tokens[1:]).replace('*', ''))
                if sace_protein:
                    cgd_infos['sace_ortholog'] = sace_protein 
                if pdb_links:
                    cgd_infos['3Dstructure'] = pdb_links
                cgd_infos['orthologs_in_non_CGD_species'] = non_cgd_orthologs
            elif text == "Best hits in Candida species" or text == "Orthologous genes in Candida species":
                for i in td.next_sibling.next_sibling.find_all('i'):
                    species = mongodatabases.get(i.get_text().strip())
                    if species:
                        ortholog_locus_tag = i.find_next('a', text=True).find_next(text=True).split('/')[0]# locus_tag of the ortholog
                        print "New Ortholog", ortholog_locus_tag, 'from entry', locus_tag
                        _id = parse_cgd_entry(client, client[species], ortholog_locus_tag, fh)# parsing of the Web page of the ortholog
                        candida_orthologs.append("%s@%s"%(_id, species))
                cgd_infos['orthologs_in_candida_species'] = candida_orthologs
            elif text == "Molecular Function" or text == "Biological Process" or text == "Cellular Component":
                for a in td.find_next('tr').find_all('a', href=True, text=True):
                    http_link = a.get('href')
                    if re.search('goid=\d+$', http_link):
                        go_definition = a.get_text().strip()
                        if 'unknown' in a.next_sibling.strip():
                            go_definition = a.next_sibling.strip().split(' (')[0]
                        go_terms ["GO:%07u"%int(http_link.split('goid=')[1])] = "%s:%s"%(go_types[text], go_definition)
                cgd_infos['go_annotations'] = go_terms
            elif text == "External Links":
                for a in td.next_sibling.next_sibling.find_all('a', href=True):
                    http_link = a.get('href')
                    if re.search('url=', http_link): # case where the http link is a redirection from the candidagenome.org to another url
                        http_link = http_link.split('url=')[1]
                    external_links.append(http_link)
                cgd_infos["external_links"] = external_links

        genomes = []
        for genome in db['genomes'].find(no_cursor_timeout = True):# timeout = False if another version of pymongo is used
            dna_molecule = DNA(name = genome['name'], sequence = genome['sequence'])
            dna_molecule.id = genome['_id']
            genomes.append(dna_molecule)

        seq_info_text = soup.find(id="alt_seq_info_1").get_text("|")
        cgd_infos['genomicStrand'] = '-' if 'Crick' in seq_info_text else '+'
        seq_info_list = seq_info_text.split("|")
        index = seq_info_list.index('GBrowse')
        cgd_infos['genomeName'] = seq_info_list[index-2].split(':')[0]

        if soup.find(id="LSPlocation"):
            features = {}
            for td in soup.find(id="LSPlocation").find_all("td", text=True):
                if td.find('a', text=True):
                    feature_type = td.find('a', text=True).find_next(text=True)
                    start = feature_type.find_next('td', text=True).find_next('td', text=True).find_next('td', text=True).find_next('td', text=True).find_next(text=True)
                    end = start.find_next('td', text=True).find_next('td', text=True).find_next(text=True)
                    coordinates = features.get(feature_type, [])
                    coordinates.append([int(start.replace(',','')), int(end.replace(',',''))] if cgd_infos['genomicStrand'] == '+' else [int(end.replace(',','')), int(start.replace(',',''))])
                    if cgd_infos['genomicStrand'] == '-':
                        coordinates = sorted(coordinates, key=lambda l: l[0])
                    features[feature_type] = coordinates
            cgd_infos['genomicPositions'] = features
        else:
            start = seq_info_list[index-2].split(':')[1].split('to')[0].strip() if cgd_infos['genomicStrand'] == '+' else seq_info_list[index-2].split(':')[1].split('to')[1].strip()
            end = seq_info_list[index-2].split(':')[1].split('to')[1].strip() if cgd_infos['genomicStrand'] == '+' else seq_info_list[index-2].split(':')[1].split('to')[0].strip()
            cgd_infos['genomicPositions'] = {cgd_infos['feature_type'].capitalize(): [[int(start), int(end)]]}
        
        for genome in genomes:
            if genome.name.split()[0] == cgd_infos['genomeName']:
                cgd_infos['genome'] = '%s@genomes'%genome.id
                for key in ['CDS', 'Noncoding_exon']+[cgd_infos['feature_type'].capitalize()]:#cgd_infos['feature_type'].capitalize() can be 'Long_terminal_repeat' or 'Repeat_region' or 'Centromere'
                    if cgd_infos['genomicPositions'].has_key(key):
                        location = cgd_infos['genomicPositions'].get(key)
                        if len(location) == 1:#gene without intron
                            cgd_infos['sequence'] = genome[location[0][0]-1:location[0][1]] if cgd_infos['genomicStrand'] == '+' else genome.get_complement()[location[0][0]-1:location[0][1]][::-1]    
                        else:#gene with intron(s)
                            spliced_sequence = ''
                            for exon in location:
                                spliced_sequence += genome[exon[0]-1:exon[1]]
                            cgd_infos['sequence'] = spliced_sequence if cgd_infos['genomicStrand'] == '+' else DNA(name=cgd_infos['locus_tag'], sequence=spliced_sequence).get_complement()[::-1]
                        break
                break

        ### creation of the MongoDB document ###
        db['annotations'].update({'_id': cgd_infos['_id']}, cgd_infos, False)
        return cgd_infos['_id']

    except (socket.timeout, socket.error), e:#OR just socket.timeout, e:
        print "With the locus %s, there was a socket error: %r"%(locus_tag, e)
        fh.write("%s\n"%locus_tag)
        annotation = db['annotations'].find_one({'locus_tag': locus_tag})
        if annotation:
            db['annotations'].remove({'locus_tag': locus_tag})
        os.execv(__file__, import_cgd(fh, locus_range), sys.argv)

if __name__ == '__main__':
    locus_range = None

    if "-locus" in sys.argv:
        locus_range = sys.argv[sys.argv.index("-locus")+1].split(',')

    dump_data()
    create_databases()
    fh = open("recovered_locus.csv", 'w')
    import_cgd(fh, locus_range)
    fh.close()
