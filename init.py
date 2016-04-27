#!/usr/bin/env python

from flask import Flask,request,render_template, jsonify
from pymongo import MongoClient
from json import dumps
import sys

sys.path.append('scripts')
import crispr_tool

app = Flask(__name__)

# Contig names for Candida glabrata CBS138
contig_names_dict = {'A': 'NC_005967',
                     'B': 'NC_005968',
                     'C': 'NC_006026',
                     'D': 'NC_006027',
                     'E': 'NC_006028',
                     'F': 'NC_006029',
                     'G': 'NC_006030',
                     'H': 'NC_006031',
                     'I': 'NC_006032',
                     'J': 'NC_006033',
                     'K': 'NC_006034',
                     'L': 'NC_006035',
                     'M': 'NC_006036',
                     'MT': 'NC_004691'}

contig_names_dict_inv = {v: k for k, v in contig_names_dict.items()}

gene_names = []
client = MongoClient('localhost', 27017)
db = client['Candida_glabrata_CBS_138']
for annotation in db['annotations'].find():
    if "standard_name" in annotation:
        gene_names.append(annotation['standard_name'])
    if "alias" in annotation:
        for alias in annotation['alias']:
            gene_names.append(alias)
    if "locus_tag" in annotation:
        gene_names.append(annotation['locus_tag'])

# We remove the duplicates
gene_names = list(set(gene_names))

@app.route('/', methods=['GET', 'POST'])
def init():
    return render_template('index.html')

@app.route('/index.html', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route('/crispr.html', methods=['GET', 'POST'])
def crispr():
    return render_template('crispr.html', data=dumps(gene_names))

@app.route('/about.html', methods=['GET', 'POST'])
def about():
    return render_template('about.html')

@app.route('/crispr_form')
def crisprform():

    result = None
    
    region_type = request.args.get('region_type', None, type=str)
    print region_type
    search_mode = 'b'
    if request.args.get('strand', None, type=str) == '+':
       search_mode = 's'
    elif request.args.get('strand', None, type=str) == '-':
       search_mode = 'a'

    pam_mode = 0
    if request.args.get('pam', None, type=str) == 'NAG':
        pam_mode = 1
    if request.args.get('pam', None, type=str) == 'Any':
        pam_mode = 2

    gc_min = 35.0
    if request.args.get('gc_min', None, type=str):
        gc_min = request.args.get('gc_min', None, type=float)

    gc_max = 75.0
    if request.args.get('gc_max', None, type=str):
        gc_max = request.args.get('gc_max', None, type=float)

    if region_type:

        # Select by gene name
        if region_type == 'collapseFour':
            gene_name = request.args.get('gene_name', None, type=str)

            result = []
            locus_tag = None
            if gene_name:
                for annotation in db['annotations'].find():
                    if "standard_name" in annotation:
                        if gene_name == annotation['standard_name']:
                            locus_tag = annotation['locus_tag']
                    if "alias" in annotation:
                        for alias in annotation['alias']:
                            if gene_name == alias:
                                locus_tag = annotation['locus_tag']
                    if "locus_tag" in annotation:
                        if gene_name == annotation['locus_tag']:
                            locus_tag = annotation['locus_tag']

                gene_in_db = False
                for grna in db['grnas'].find():
                    if grna['gene_name'] == locus_tag:
                        gene_in_db = True
                        gcContent = grna['gcContent']
                        if gcContent >= gc_min and gcContent <= gc_max:
                            pamSequence = grna['pamSequence']
                            if (pam_mode == 0 and pamSequence[1:] == 'GG') or (pam_mode == 1 and pamSequence[1:] == 'AG') or (pam_mode == 2 and pamSequence[2:] == 'G'): 
                                genomicStrand = grna['genomicStrand']
                                if (search_mode == 's' and genomicStrand == '+') or (search_mode == 'a' and genomicStrand == '-') or search_mode == 'b':
                                    result.append(grna)
                if gene_in_db == False:
                    result = crispr_tool.find(db_name='Candida_glabrata_CBS_138', genome_name=None, start=None, end=None, pam_mode=pam_mode, search_mode=search_mode, gc_min=gc_min, gc_max=gc_max, mism=5, guides_file=None, host='localhost', port=27017, gene_name=gene_name, sequence=None)

        # Select by coordinates
        elif region_type == 'collapseFive':
            chromosome = request.args.get('chromosome', None, type=str)
            genomic_start = request.args.get('genomic_start', None, type=int)
            genomic_end = request.args.get('genomic_end', None, type=int)
            if chromosome and genomic_start and genomic_end:
                print "pam :", pam_mode
                result = crispr_tool.find(db_name='Candida_glabrata_CBS_138', genome_name=chromosome, start=genomic_start, end=genomic_end, pam_mode=pam_mode, search_mode=search_mode, gc_min=gc_min, gc_max=gc_max, mism=5, guides_file=None, host='localhost', port=27017, gene_name=None, sequence=None)
        # Select by sequence
        elif region_type == 'collapseSix':
            sequence = request.args.get('sequence', None, type=str)
            if sequence:
                result = crispr_tool.find(db_name='Candida_glabrata_CBS_138', genome_name=None, start=None, end=None, pam_mode=pam_mode, search_mode=search_mode, gc_min=gc_min, gc_max=gc_max, mism=5, guides_file=None, host='localhost', port=27017, gene_name=None, sequence=sequence)

    return dumps(result)

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

if __name__ == '__main__':
    app.run() #if app.run() => localhost by default