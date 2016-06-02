#!/usr/bin/env python

import os, re, shutil, commands
from pyrna.parsers import to_fasta
from pyrna.features import DNA
from pyrna.computations import Tool, Blast
from pyrna import utils
from pymongo import MongoClient
from pandas import DataFrame
import purge_dir, purge_dir_blastdb
from Bio import AlignIO

"""
COMMAND-LINE
------------
You can run the program with this command-line:
./search_orthologs.py
"""

class Hmmer(Tool):
    def __init__(self, seqdb, cache_dir="/tmp"):
        Tool.__init__(self, cache_dir = cache_dir)
        self.seqdb = seqdb

    def hmmbuild(self, multiple_alignment, query_name):
        """
        Constructs a profile HMM from a multiple sequence alignement (DNA or proteins).
        Format CLUSTAL W is read.
        """
        path = self.cache_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)
        with open("%s/input.aln"%path, 'w+b') as clustal_file:
            clustal_file.write(multiple_alignment)
        
        self.hmm_profile = path+"/profile.hmm"
        commands.getoutput("cd %s ; hmmbuild -n %s %s %s"%(self.cache_dir, query_name, self.hmm_profile, clustal_file.name))

    def hmmsearch(self):
        """
        hmmsearch searches profile(s) against a sequence database.

        hmmsearch [options] <hmmfile> <seqdb>
        <seqdb> BLAST database 
        """
        return self.parse_hmmsearch_output(commands.getoutput("cd %s ; hmmsearch %s %s"%(self.cache_dir, self.hmm_profile, self.seqdb)))

    def parse_hmmsearch_output(self, output):
        """
        Parses the hmmsearch output.

        Parameters:
        ---------
        - output: the hmmsearch output content as a String 

        Returns:
        --------
        A pandas DataFrame describing all the hmmsearch hits. The index stores hit ids. The columns are:
        - query
        - target
        - e_value
        - score (bit score)
        """
        hits=[]
        query_name = None
        target_name = None
        evalue = None
        tag = False
        lines = output.split('\n')
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line:
                if line.startswith("Query:"):
                    query_name = line.split()[1]
                elif line.startswith("E-value") and not tag:
                    tag = True
                elif not line.startswith("-") and not line.startswith("Domain") and not line.startswith("[No") and tag:
                    tokens = line.split()
                    hits.append({
                        "query": query_name,
                        "e_value": float(tokens[0]),
                        "score": float(tokens[1]),
                        "target": tokens[8]
                        })
                elif line.startswith("[No"):
                    hits.append({
                        "query": query_name,
                        "e_value": '',
                        "score": '',
                        "target": ''
                        })
                    break 
                elif (line.startswith("-") or line.startswith("Domain")) and hits:
                    break
            i += 1
        return DataFrame(hits)

class Muscle(Tool):
    def __init__(self, fasta_sequences, cache_dir="/tmp"):
        Tool.__init__(self, cache_dir = cache_dir)
        
        path = cache_dir+"/"+utils.generate_random_name(7)
        os.mkdir(path)
        with open("%s/input.fasta"%path, 'w+b') as fasta_file:
            fasta_file.write(fasta_sequences)

        self.fasta_sequences = fasta_file.name

    def align(self):
        """
        Aligns multiple sequences using the MUSCLE algorithm. 
        """
        return self.remove_header(commands.getoutput("cd %s ; muscle -in %s"%(self.cache_dir, self.fasta_sequences)))

    def realign(self, existing_msa, new_seqs=None):
        """
        Aligns a multiple sequence alignment with new sequence(s) using the MUSCLE algorithm.
        Command line: muscle -profile -in1 existing_aln.afa -in2 new_seqs.afa 
        The output alignment will be at the CLUSTAL W format and will contain a CONSENSUS line. If the output format is FASTA, no consensus line ! 
        """
        tmp_dir = os.path.dirname(self.fasta_sequences)
        with open("%s/existing_msa.fasta"%tmp_dir, 'w+b') as existing_msa_file:
            existing_msa_file.write(existing_msa)

        if new_seqs:
            with open("%s/new_seqs.fasta"%tmp_dir, 'w+b') as new_seqs_msa_file:
                new_seqs_msa_file.write(new_seqs)
            return self.remove_header(commands.getoutput("cd %s ; muscle -profile -clwstrict -in1 %s -in2 %s"%(self.cache_dir, existing_msa_file.name, new_seqs_msa_file.name)))
        else:
            return self.remove_header(commands.getoutput("cd %s ; muscle -profile -clwstrict -in1 %s -in2 %s"%(self.cache_dir, existing_msa_file.name, self.fasta_sequences)))

    def remove_header(self, output):
        """
        Removes the header of a Muscle output at the format CLUSTAL W.
        """
        new_output = []
        tag = False
        lines = output.split('\n')
        for line in lines:
            if line.startswith(">") or line.startswith("CLUSTAL W (1.81) multiple sequence alignment") or line.startswith("MUSCLE (3.8) multiple sequence alignment"):
                tag = True
            if tag:   
                new_output.append(line)
        return '\n'.join(new_output)    

def convert_msa_format(clustal_multiple_alignment, input_format, output_format):
    """
    Function that converts the format of a multiple sequence alignment (MSA).
    """
    path = "/tmp/"+utils.generate_random_name(7)
    os.mkdir(path)
    with open("%s/msa.aln"%path, 'w+b') as clustal_file:
        clustal_file.write(clustal_multiple_alignment)
    input_infh = open("%s/msa.aln"%path, 'rU')
    output_outfh = open("%s/msa.fasta"%path, 'w')
    alignment = AlignIO.parse(input_infh, input_format)
    AlignIO.write(alignment, output_outfh, output_format)
    output_outfh.close()
    input_infh.close()
    with open("%s/msa.fasta"%path, 'r') as output_infh:
        fasta_multiple_alignment = output_infh.read()
        return fasta_multiple_alignment

def parse_protein_clustalw(clustalw_data):
    """
    Parses protein data at the format CLUSTAL W.

    Parameters:
    ---------
     - clustalw_data: the CLUSTAL W protein data as a String

    Returns:
    ------
    A list of PyRNA DNA objects
    """
    alignedSequences = {}
    lines = clustalw_data.strip().split('\n')
    for line in lines:
        line = line.strip()
        if len(line) and not line.startswith('CLUSTAL W'):
            if not line.startswith('*') and not line.startswith('.') and not line.startswith(':'):
                tokens = line.split()
                alignedSequences[tokens[0]] = alignedSequences.get(tokens[0], "") + tokens[1]
            else:# consensus sequence
                pass

    proteins = []

    for key in alignedSequences:
        prot = DNA(name=key.split('/')[0], sequence=alignedSequences[key])
        proteins.append(prot)

    return proteins

def purge_dir(directory, option=None):
    """
    Removes in a given directory (for example "/tmp") the sub-directories created by the script search_orthologs.py

    Parameters:
    ---------
    - directory
    - option: None or "blastdb"
    With the option "blastdb", we remove sub-directories created by the command 'formatdb' (formats a database that will be used by the tool BLASTP).
    """
    directories = os.listdir(directory)

    tag = False
    for d in directories:
        if os.path.isdir(directory+'/'+d): # test if d is a directory
            if re.match('^\w+$', d) and os.access(directory+'/'+d, os.R_OK):
                files = os.listdir(directory+'/'+d)
                for f in files:
                    if option == 'blastdb':
                        if f.endswith('.pin'): # file for a database used by BLASTP
                            tag = True    
                    else:
                        if f.endswith('.hmm') or f.startswith('existing_msa.') or f.startswith('msa.'):# file for a database used by Hmmer() or Muscle() or convert_msa_format()
                            tag = True
            if tag:
                shutil.rmtree(directory+'/'+d)
                tag = False

def hmm():
    """
    Searches orthologs of the Candida species in the Nakaseomyces species.
    From alignments of multiple protein sequences of the Candida species, building of HMM profils that will allow to find 1 or several orthologs in the Nakaseomyces species using the tool "hmmsearch".
    
    Step 1: building of formatted BLAST databases from protein sequences of the Nakaseomyces species.

    Step 2: for each multiple sequence alignment of Candida proteins, building of an HMM profil. Then, searching of Nakaseomyces orthologs using the tool "hmmsearch".
        Utilization of the bit score to retrieve the orthologs having a score >= 70 % of the best score.

    Step 3: realignment with Muscle (tool used to do the original alignment on the CGD Website): alignment of the Nakaseomyces orthologs with all the Candida orthologs

    Step 4: updating data in the Mongo databases of the Candida species: 
        table: 'annotations'
        field: 'orthologs_in_candida_species' => cgd_orthologs + nakaseo_orthologs

    Step 5: updating data in the Mongo databases of the Nakaseomyces species:
        table: 'annotations'
        field: 'alignment' => alignment_id
        field: 'orthologs_in_candida_species' => cgd_orthologs + nakaseo_orthologs
    
    Step 6: updating data in the Mongo database named 'comparative_genomics':
        table: 'proteins'
        field: 'alignment' => new alignment
    """
    client = MongoClient()

    cgd_dbs = ['Candida_glabrata_CBS_138', 'Candida_albicans_SC5314', 'Candida_dubliniensis_CD36', 'Candida_parapsilosis_CDC317']
    cgd_dict ={}
    for species in cgd_dbs:
        cgd_molecules = []
        for annotation in client[species]['annotations'].find({'translation':{'$exists':True}}, timeout = False):
            dna = DNA(name=annotation['locus_tag'], sequence=annotation['translation'])
            dna.id = annotation['_id']
            dna.source = annotation['source']
            cgd_molecules.append(dna)
        cgd_dict[species] = cgd_molecules

    ### CREATION OF FORMATTED BLAST DATABASES FROM PROTEINS OF THE NAKASEOMYCES SPECIES ###
    nakaseo_dbs = ['Nakaseomyces_bracarensis_CBS_10154', 'Nakaseomyces_castellii_CBS_4332', 'Nakaseomyces_nivariensis_CBS_9983', 'Nakaseomyces_bacillisporus_CBS_7720', 'Nakaseomyces_delphensis_CBS_2170']
    nakaseo_dict = {}
    for species in nakaseo_dbs:
        nakaseo_molecules = []
        for annotation in client[species]['annotations'].find(timeout = False):
            dna = DNA(name=annotation['locus_tag'], sequence=annotation['translation'])
            dna.id = annotation['_id']
            dna.source = annotation['source'] # for example: "db:gryc:Candida_bracarensis_CBS_10154"
            nakaseo_molecules.append(dna)
        blast = Blast(target_molecules=nakaseo_molecules)
        blast.format_db(is_nucleotide=False) # use formatdb of RNA_algo/blast-2.2.26 (not Tools/ncbi-blast-2.2.29+)
        nakaseo_dict[species] = [blast.formatted_db, nakaseo_molecules]

    ### FOR EACH PROTEIN ALIGNMENT, SEACHING OF ORTHOLOGS IN THE NAKASEOMYCES SPECIES AND REALIGNMENT ###
    aln_total = client['comparative_genomics']['proteins'].find().count()
    print "Total number of aligned loci: %i"%aln_total 
    aln_counter = 0
    hmm_results = {} 
    for protein in client['comparative_genomics']['proteins'].find(timeout = False):
        aln_counter += 1
        locus_tag = protein['locus_tag']
        alignment_id = protein['_id']
        print "locus number %i named %s in progress (%.2f %% of total aligned loci)"%(aln_counter, locus_tag, (aln_counter/float(aln_total))*100)
        multiple_alignment = protein['alignment'] # multiple protein sequence alignment at CLUSTAL W format
        clustalw_molecules = parse_protein_clustalw(multiple_alignment)
        
        fasta_alignment = convert_msa_format(multiple_alignment, "clustal", "fasta")

        nakaseo_to_align = []
        nakaseo_orthologs = []
        for species in nakaseo_dbs:
            hmmer = Hmmer(seqdb=nakaseo_dict[species][0]) # seqdb = BLAST protein database

            hmmer.hmmbuild(multiple_alignment, locus_tag)
            df = hmmer.hmmsearch()
            if not df.empty:
                i = 0
                best_score = None
                for row in df.iterrows():
                    i += 1
                    data_hmmer = row[1]
                    if data_hmmer['target']:
                        if i == 1:
                            best_score = data_hmmer['score']
                            percent = 100
                        elif i > 1:
                            percent = (data_hmmer['score'] / best_score) * 100
                        if percent >= 70:
                            for nakaseo_molecule in nakaseo_dict[species][1]:
                                if nakaseo_molecule.name == data_hmmer['target']:
                                    if nakaseo_molecule.sequence.endswith('*'):
                                        nakaseo_molecule.sequence = nakaseo_molecule.sequence[:-1]
                                    nakaseo_to_align.append(nakaseo_molecule)
                                    nakaseo_orthologs.append("%s@%s"%(nakaseo_molecule.id, nakaseo_molecule.source.split(':')[-1]))# for example, nakaseo_molecule.source = "db:gryc:Candida_bracarensis"
                                    break

        ### SOME ORTHOLOGS OF THE CANDIDA SPECIES ARE NOT PRESENT IN THE ORIGINAL ALIGNMENT ###
        cgd_to_align = []
        cgd_orthologs = []
        for cgd_db in cgd_dbs: 
            for annotation in client[cgd_db]['annotations'].find({'alignment': alignment_id}, timeout = False): # several annotations (different locus_tags) can have the same alignment
                cgd_orthologs.append("%s@%s"%(annotation['_id'], cgd_db))
                for cgd_molecule in cgd_dict[cgd_db]:
                    if cgd_molecule.name == annotation['locus_tag']:
                        tag = True
                        for clustalw_molecule in clustalw_molecules:
                            if cgd_molecule.name == clustalw_molecule.name:
                                tag = False
                        if tag:
                            cgd_to_align.append(cgd_molecule)
                        break

        molecules_to_align = cgd_to_align+nakaseo_to_align 
        
        ### REALIGNMENT WITH ALL ORTHOLOGOUS SEQUENCES ###
        if molecules_to_align:
            muscle = Muscle(fasta_sequences=to_fasta(molecules=molecules_to_align))
            if len(molecules_to_align) == 1:
                output_msa = muscle.realign(existing_msa=fasta_alignment)
            else:
                input_msa = muscle.align()
                output_msa = muscle.realign(existing_msa=fasta_alignment, new_seqs=input_msa)

            ### USEFUL FOR THE FOLLOWING STEP OF MSA UPDATING ###
            
            hmm_results[alignment_id] = output_msa
            
            ### DATA UPDATING IN THE CANDIDA SPECIES MONGO DATABASES ###
            for cgd_db in cgd_dbs:
                for annotation in client[cgd_db]['annotations'].find({'alignment': alignment_id}, timeout = False):
                    cgd_orthologs_copy = cgd_orthologs[:]
                    for ortholog in cgd_orthologs:
                        if annotation['_id'] in ortholog:
                            cgd_orthologs_copy.remove(ortholog) # removes the current CGD ortholog in the list   
                    client[cgd_db]['annotations'].update({'_id': annotation['_id']},{'$set':{'orthologs_in_candida_species': cgd_orthologs_copy+nakaseo_orthologs}}, False)     
            
            ### DATA UPDATING IN THE NAKASEOMYCES SPECIES MONGO DATABASES ###
            for nakaseo in nakaseo_to_align:
                nakaseo_db = nakaseo.source.split(':')[-1]
                nakaseo_orthologs_copy = nakaseo_orthologs[:]
                for ortholog in nakaseo_orthologs:
                    if nakaseo.id in ortholog:
                        nakaseo_orthologs_copy.remove(ortholog) # removes the current nakaseo ortholog in the list     
                client[nakaseo_db]['annotations'].update({'_id': nakaseo.id},{'$set':{'alignment': alignment_id, 'orthologs_in_candida_species': cgd_orthologs+nakaseo_orthologs_copy}}, False)
    
        purge_dir(directory="/tmp")
    purge_dir(directory="/tmp", option="blastdb")

    ### MSA UPDATING IN THE COMPARATIVE GENOMICS MONGO DATABASE ###
    for item in hmm_results.iteritems(): 
        client['comparative_genomics']['proteins'].update({'_id': item[0]},{'$set':{'alignment': item[1]}}, False)
    
    client.disconnect()

if __name__ == '__main__':

    hmm() 
