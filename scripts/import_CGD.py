#!/usr/bin/env python

import os, sys, re, urllib, urllib2
from pyrna.parsers import parse_fasta
from pyrna.features import Molecule
from pymongo import MongoClient
from bson import ObjectId
from bs4 import BeautifulSoup

"""
COMMAND-LINE
------------
You can run the program with this command-line:
./construct_db.py
"""

def dump_data():
    """
    This function recovers from the database GRYC (http://www.igenolevures.org) all the data necessary to construct the CAGGLE DB (except for C. dubliniensis which is recovered from NCBI). It will produce one directory per species. In each directory, it will dump the data in these files:
    - sequences.fasta: genomic sequences
    - proteins.fasta: protein sequences.
    For the species C. glabrata, this script will also produce a file named locus_tag.txt. This file will be used for the next step (see function import_cgd())
    """

    dirname = os.path.abspath(os.path.dirname(__file__))+'/../data'
    if not os.path.exists(dirname):
        #we create the data dir
        print "dir not exist"
    species = ['Candida_glabrata_CBS_138'] #the species we're interested in
    for s in species:
        if not os.path.exists(dirname+'/'+s):
            #we create the species dir
            pass
        #now we recover and dump the data in this dir

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
        'class': "CDS"
        'locus_tag'
        'translation': protein sequence
        'description': gene annotation of GRYC
        'genome': "genome_id@genomes"
        'genomeName'
    """
    client = MongoClient()

    dirname = os.path.abspath(os.path.dirname(__file__))+'/../data'

    print dirname
    pattern = re.compile('^BN\d{3,3}(_\w{4,4}\ds\d\d)e.+$')

    for species in os.listdir(dirname):
        if os.path.isdir("%s/%s"%(dirname, species)):
            if "sequences.fasta" in os.listdir(dirname+'/'+species):
                with open("%s/%s/sequences.fasta"%(dirname, species)) as h:
                    content = h.read()
                    dna_molecules = parse_fasta(content, 'DNA')
                    for d in dna_molecules:
                        client[species]['genomes'].insert(
                            {
                                '_id': str(ObjectId()),
                                'name': d.name,
                                'sequence': d.sequence
                            }
                        )
            if "proteins.fasta" in os.listdir(dirname+'/'+species):
                with open("%s/%s/proteins.fasta"%(dirname, species)) as h:
                    content = h.read()
                    protein_molecules = parse_fasta(content, 'DNA')
                    for p in protein_molecules:
                        tokens = p.name.split()
                        if p.sequence.endswith("*"):
                            p.sequence = p.sequence[:-1]
                        infos = {
                                '_id': str(ObjectId()),
                                'source': "db:gryc:%s"%species,
                                'class': 'CDS',
                                'locus_tag': tokens[0],
                                'translation': p.sequence,
                                'description': ' '.join(tokens[1:])
                                }
                        match = pattern.search(tokens[0])
                        if match:
                            infos['genomeName'] = match.group(1)
                            genome = client[species]['genomes'].find_one({'name': infos['genomeName']})
                            if genome:
                                infos['genome'] = genome['_id']+'@genomes'
                        else:
                            print "Warning: locus tag %s has no genome name!"%tokens[0]
                        client[species]['annotations'].insert(infos)

    client.disconnect()

def import_cgd():
    """
    From each gene of Candida glabrata CBS138, some informations are extracted from its HTML page on the site CGD (Candida Genome Database).

    Step 1. Creation of a MongoDB with the name 'comparative_genomics' that contains a table named 'proteins' with the following fields:
        '_id': ID of the alignment
        'locus_tag': locus tag of the Candida glabrata gene
        'alignment': Clustalw multiple alignment of protein sequences

    Step 2. Creation of a table named 'annotations' for the four CGD Candida species with the following fields:
        '_id'
        'locus_tag'
        'source': "db:cgd:locus_tag"
        'standard_name'
        'alias'
        'genomicsStrand'
        'genome': "genome_id@genomes"
        'genomeName'
        'genomicPositions'
        'feature_type'
        'description'
        'orthologs_in_candida_species': list of "annotation_id@species_name"
        'orthologs_in_non_CGD_species': list of HTTP links
        'go_annotations'
        'external_links'
        'translation'
        'alignment': ID of the multiple alignment of protein sequences
    """

    client = MongoClient()

    db = client['Candida_glabrata_CBS_138']

    dirname = os.path.abspath(os.path.dirname(__file__))
    with open("%s/Candida_glabrata_CBS_138/locus_tags.txt"%dirname) as h:
        for locus_tag in h.readlines():
            parse_cgd_entry(client, db, locus_tag.strip())# parse_cgd_entry() is a recursive function

    client.disconnect()

def parse_cgd_entry(client, db, locus_tag, protein_alignment_id = None):

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
    soup = BeautifulSoup(str(urllib.urlopen(url).read()))
    for a in soup.find_all(id="tabtext", href=True):
        if a.next_sibling == ' Protein':
            protein_soup = BeautifulSoup(str(urllib.urlopen(a.get('href')).read()))
            for td in protein_soup.find_all("td", {"bgcolor":"#FFCC33"}):
                if td.get_text().strip() == 'in FASTA format)':
                    params = urllib.urlencode({'downloadFasta': '1', 'feature_name': td.find("form").find_next('input').find_next('input').get('value')})
                    response = urllib2.urlopen("http://www.candidagenome.org%s"%td.find("form").get('action'), params)
                    fasta_protein = response.read()
                    cgd_infos['translation'] = parse_fasta(fasta_protein, 'DNA')[0].sequence.replace('*', '')
        elif a.next_sibling == ' Homologs':
            if locus_tag.startswith('CAGL') or locus_tag.startswith('Novel_') or locus_tag.startswith('Cagl'):
                homologs_soup = BeautifulSoup(str(urllib.urlopen(a.get('href')).read()))
                for a in homologs_soup.find_all("a", {"target":"infowin"}):
                    if a.get_text().strip() == 'Protein alignment (ClustalW format)' and not client['comparative_genomics']['proteins'].find_one({'locus_tag':locus_tag}):
                        protein_alignment_id = str(ObjectId())
                        client['comparative_genomics']['proteins'].insert( {
                            '_id': protein_alignment_id,
                            'locus_tag': locus_tag,
                            'alignment': str(urllib.urlopen(a.get('href')).read())
                        })
    for td in soup.find_all("td", {"bgcolor":"#FFCC33", 'valign':'top'}):
        text = td.get_text().strip()
        if text == "Standard Name":
            cgd_infos['standard_name'] = td.next_sibling.next_sibling.get_text("|").split("|")[0]
        elif text == "Systematic Name, Reference Strain":
            """
            To eliminate redundancy with the locus tags in the Candida albicans MongoDB
            For example: Systematic Name 'C1_01530C_A' = Assembly 19/21 Identifier 'orf19.3341'
            """
            if locus_tag.startswith('orf') and db.name == 'Candida_albicans_SC5314':
                locus_tag = td.next_sibling.next_sibling.get_text().split()[0]
                cgd_infos['locus_tag'] = locus_tag
                cgd_infos['source'] = 'db:cgd:%s'%locus_tag
                cgd_annotation = db['annotations'].find_one({'locus_tag': locus_tag})
                if cgd_annotation:
                    print locus_tag, "already stored"
                    return cgd_annotation['_id']
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
            for i in td.next_sibling.next_sibling.find_all('i'):
                species = i.get_text().strip()
                http_link = i.find_next('a', href=True).get('href')
                non_cgd_orthologs.append("%s:%s"%(species, http_link))
            cgd_infos['orthologs_in_non_CGD_species'] = non_cgd_orthologs
        elif text == "Best hits in Candida species" or text == "Orthologous genes in Candida species":
            for i in td.next_sibling.next_sibling.find_all('i'):
                species = mongodatabases.get(i.get_text().strip())
                if species:
                    ortholog_locus_tag = i.find_next('a', text=True).find_next(text=True).split('/')[0]# locus_tag of the ortholog
                    print "New Ortholog", ortholog_locus_tag, 'from entry', locus_tag
                    _id = parse_cgd_entry(client, client[species], ortholog_locus_tag, protein_alignment_id)# parsing of the Web page of the ortholog
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

    candida_chromosomes = {
                    'ChrA_C_glabrata_CBS138': 'NC_005967',
                    'ChrB_C_glabrata_CBS138': 'NC_005968',
                    'ChrC_C_glabrata_CBS138': 'NC_006026',
                    'ChrD_C_glabrata_CBS138': 'NC_006027',
                    'ChrE_C_glabrata_CBS138': 'NC_006028',
                    'ChrF_C_glabrata_CBS138': 'NC_006029',
                    'ChrG_C_glabrata_CBS138': 'NC_006030',
                    'ChrH_C_glabrata_CBS138': 'NC_006031',
                    'ChrI_C_glabrata_CBS138': 'NC_006032',
                    'ChrJ_C_glabrata_CBS138': 'NC_006033',
                    'ChrK_C_glabrata_CBS138': 'NC_006034',
                    'ChrL_C_glabrata_CBS138': 'NC_006035',
                    'ChrM_C_glabrata_CBS138': 'NC_006036',
                    'mito_C_glabrata_CBS138': 'NC_004691'
                    }

    genomes = []
    for genome in db['genomes'].find(timeout = False):
        dna_molecule = Molecule(name = genome['name'])
        dna_molecule.id = genome['_id']
        dna_molecule.sequence = genome['sequence']
        genomes.append(dna_molecule)

    seq_info_text = soup.find(id="alt_seq_info_1").get_text("|")
    cgd_infos['genomicStrand'] = '-' if 'Crick' in seq_info_text else '+'
    seq_info_list = seq_info_text.split("|")
    index = seq_info_list.index('GBrowse')
    cgd_infos['genomeName'] = candida_chromosomes[seq_info_list[index-2].split(':')[0]] if db.name == 'Candida_glabrata_CBS_138' else seq_info_list[index-2].split(':')[0]

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

    for genome in genomes:
        if genome.name.split()[0] == cgd_infos['genomeName']:
            cgd_infos['genome'] = '%s@genomes'%genome.id
            break

    ### creation of the MongoDB document ###
    if protein_alignment_id:
        cgd_infos['alignment'] = protein_alignment_id
    db['annotations'].update({'_id': cgd_infos['_id']}, cgd_infos, False)

    return cgd_infos['_id']

if __name__ == '__main__':
    dump_data()
    create_databases()
    import_cgd()
