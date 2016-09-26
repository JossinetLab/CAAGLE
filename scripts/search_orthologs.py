#!/usr/bin/env python

import sys, os, re, shutil, commands, scipy.stats
from pyrna.parsers import to_fasta
from pyrna.features import DNA
from pyrna.computations import Tool, Blast
from pyrna import utils
from pymongo import MongoClient
from pandas import DataFrame
from Bio import AlignIO
from bson import ObjectId

"""
COMMAND-LINE
------------
You can run the program with this command-line:
./search_orthologs.py
"""

class MyBlast(Blast):
    def __init__(self, target_molecules, cache_dir):# cache_dir="/tmp" in the original class Blast
        Blast.__init__(self, target_molecules, cache_dir)

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
        # print output
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

    def align(self, output_format=None):
        """
        Aligns multiple sequences using the MUSCLE algorithm. 
        """
        format = ''#FASTA format
        if output_format == "clustalw":
            format = ' -clwstrict'

        return self.remove_header(commands.getoutput("cd %s ; muscle -in %s%s"%(self.cache_dir, self.fasta_sequences, format)))#IF HEADER

    def realign(self, existing_msa, new_seqs=None, output_format=None):
        """
        Aligns a multiple sequence alignment with new sequence(s) using the MUSCLE algorithm.
        Command line: muscle -profile -in1 existing_aln.afa -in2 new_seqs.afa 
        The output alignment will be at the FASTA (by default) or CLUSTAL W format and will contain a CONSENSUS line. If the output format is FASTA, no consensus line ! 
        """
        tmp_dir = os.path.dirname(self.fasta_sequences)
        with open("%s/existing_msa.fasta"%tmp_dir, 'w+b') as existing_msa_file:
            existing_msa_file.write(existing_msa)

        format = ''#FASTA format
        if output_format == "clustalw":
            format = ' -clwstrict '

        if new_seqs:
            with open("%s/new_seqs.fasta"%tmp_dir, 'w+b') as new_seqs_msa_file:
                new_seqs_msa_file.write(new_seqs)
            return self.remove_header(commands.getoutput("cd %s ; muscle -profile%s -in1 %s -in2 %s"%(self.cache_dir, format, existing_msa_file.name, new_seqs_msa_file.name)))
        else:
            return self.remove_header(commands.getoutput("cd %s ; muscle -profile%s -in1 %s -in2 %s"%(self.cache_dir, format, existing_msa_file.name, self.fasta_sequences)))

    def get_spscore(self):
        spscore_output = commands.getoutput("cd %s ; muscle -spscore %s"%(self.cache_dir, self.fasta_sequences))
        spscore = None
        lines = spscore_output.split('\n')
        for line in lines:
            if line.startswith('File='):#'File=alignment.fsa;SP=18.97'
                spscore = float(line.split('SP=')[1].strip())
        return spscore

    def remove_header(self, msa_output):
        """
        Removes the header of a Muscle output.
        """
        msa_wo_header = []
        tag = False
        lines = msa_output.split('\n')
        for line in lines:
            if line.startswith(">") or line.startswith("CLUSTAL W (1.81) multiple sequence alignment") or line.startswith("MUSCLE (3.8) multiple sequence alignment"):
                tag = True
            if tag:   
                msa_wo_header.append(line)
        return '\n'.join(msa_wo_header)    

def convert_msa_format(input_msa, input_format, output_format):
    """
    Function that converts the format of a multiple sequence alignment (MSA).
    """
    path = "/tmp/"+utils.generate_random_name(7)
    os.mkdir(path)
    with open("%s/msa.aln"%path, 'w+b') as msa_file:
        msa_file.write(input_msa)
    input_infh = open("%s/msa.aln"%path, 'rU')
    output_outfh = open("%s/msa.fasta"%path, 'w')
    alignment = AlignIO.parse(input_infh, input_format)
    AlignIO.write(alignment, output_outfh, output_format)
    output_outfh.close()
    input_infh.close()
    with open("%s/msa.fasta"%path, 'r') as output_infh:
        output_msa = output_infh.read()
        return output_msa

def purge_dir(directory, option=None):
    """
    Removes in a given directory the sub-directories created by the script search_orthologs.py

    Parameters:
    ---------
    - directory
    - option: None or "blastdb"
    With the option "blastdb", we remove sub-directories created by the command 'formatdb' (formats a database that will be used by the tool BLASTP).
    """
    directories = os.listdir(directory)

    tag = False
    for d in directories:
        if os.path.isdir(directory+'/'+d): # we test if d is a directory
            if re.match('^\w+$', d) and os.access(directory+'/'+d, os.R_OK):
                files = os.listdir(directory+'/'+d)
                for f in files:
                    if option == 'blastdb':
                        if f.endswith('.pin'): # file for a database used by BLASTP
                            tag = True    
                    else:
                        if not d.endswith('_db'):
                            if f.endswith('.hmm') or f.endswith('.fasta') or f.endswith('.aln'):# file for a database used by the functions Hmmer(), Muscle() and convert_msa_format()
                                tag = True
            if tag:
                shutil.rmtree(directory+'/'+d)
                tag = False

def hmm(hmm_threshold):
    """
    Arguments:
    --------- 
    -hmm_threshold: Integer (by default: 100 ; the best score) threshold percent / best bit score of hmmsearch

    Description:
    -----------
    This function aligns orthologs of the Candida and/or S. cerevisiae species and searches orthologs in the Nakaseomyces species.
    
    From each gene of C. glabrata and its Candida and/or S. cerevisiae orthologs, it creates an alignment of multiple protein sequences and builds a HMM profil that will allow to find orthologs in the Nakaseomyces species using the tool "hmmsearch".
    
    Step 1: alignment with the Muscle algorithm of the protein sequences of a C. glabrata gene and its Candida and/or S. cerevisiae orthologous genes
    Step 2: building formatted BLAST databases from protein sequences of the Nakaseomyces species.
    Step 3: for each multiple sequence alignment of Candida proteins, building an HMM profil. Then, searching Nakaseomyces orthologs using the tool "hmmsearch".
        Utilization of the bit score to retrieve the orthologs having a score >= n percent of the best score (see the argument -hmm_threshold).
    Step 4: if nakaseo orthologs exist, realignment with Muscle i.e. alignment of the Nakaseomyces sequences with its Candida and/or cerevisiae orthologous sequences
            otherwise, the original alignment is stored in the MongoDB
    Step 5: updating data in the Mongo databases of the Candida and Nakaseomyces species:
        table 'annotations'
            field 'alignment': list of alignment IDs ; for C. glabrata, this list contains a unique ID
            field 'orthologs_in_candida_species': list of 'annotations_id@database_name'
    Step 6: creation of a Mongo database named 'comparative_genomics':
        table 'proteins'
            field '_id': proteins ID = alignment ID
            field 'locus_tag': name of the C. glabrata locus used to do the original alignment
            field 'alignment': multiple sequence alignments = dictionary with keys 'all_species' AND 'non_pathogenic_species' AND/OR 'pathogenic_species'
            field 'spscore': SP scores of the MSA = dictionary with keys 'all_species' AND 'non_pathogenic_species' AND/OR 'pathogenic_species'
            field 'percentile_of_spscore': percentiles of SP scores =  dictionary with keys 'all_species' AND 'non_pathogenic_species' AND/OR 'pathogenic_species'

    Comments:
        Species in CGD = ['Candida_glabrata_CBS_138', 'Candida_albicans_SC5314', 'Candida_dubliniensis_CD36' and 'Candida_parapsilosis_CDC317']
        Pathogenic species = ['Candida_glabrata_CBS_138', 'Candida_albicans_SC5314', 'Candida_dubliniensis_CD36', 'Candida_parapsilosis_CDC317', 'Nakaseomyces_bracarensis_CBS_10154', 'Nakaseomyces_nivariensis_CBS_9983'] 
    """
    non_pathogens = ['Nakaseomyces_bacillisporus_CBS_7720', 'Nakaseomyces_castellii_CBS_4332', 'Nakaseomyces_delphensis_CBS_2170']#and Saccharomyces cerevisiae S288C

    client = MongoClient()

    all_spscores = []
    non_patho_spscores = []
    patho_spscores = []

    ### CREATION OF FORMATTED BLAST DATABASES FROM PROTEINS OF THE NAKASEOMYCES SPECIES ###
    nakaseo_dbs = ['Nakaseomyces_bracarensis_CBS_10154', 'Nakaseomyces_castellii_CBS_4332', 'Nakaseomyces_nivariensis_CBS_9983', 'Nakaseomyces_delphensis_CBS_2170', 'Nakaseomyces_bacillisporus_CBS_7720']
    nakaseo_dict = {}
    for species in nakaseo_dbs:
        nakaseo_molecules = []
        for annotation in client[species]['annotations'].find({'translation':{'$exists':True}}, no_cursor_timeout = True):#tRNA, ncRNA and rRNA have not translated sequence
            dna = DNA(name=annotation['locus_tag'], sequence=annotation['translation'])
            dna.id = annotation['_id']
            nakaseo_molecules.append(dna)
        blast = Blast(target_molecules=nakaseo_molecules, cache_dir="/tmp/%s_db"%species)
        blast.format_db(is_nucleotide=False)
        nakaseo_dict[species] = [blast.formatted_db, nakaseo_molecules]

    total_loci = client['Candida_glabrata_CBS_138']['annotations'].find({'alignment':{'$exists':False}}).count()
    print "Total number of genes processed: %i"%total_loci
    loci_counter = 0

    ### FROM EACH C. GLABRATA LOCUS THAT HAS NOT ALIGNMENT IN THE MONGODB ###
    for glabrata_annotation in client['Candida_glabrata_CBS_138']['annotations'].find({'alignment':{'$exists':False}}, no_cursor_timeout = True):#find({'locus_tag':{'$in':["CAGL0M05665g", "CAGL0G06006g"]}}):
        loci_counter += 1
        locus_tag = glabrata_annotation['locus_tag']
        print "Locus number %i named %s in progress (%.2f %% of total loci)"%(loci_counter, locus_tag, (loci_counter/float(total_loci))*100)
        aligned_orthologs = []
        cgd_molecules_to_align = []
        nakaseo_molecules_to_align = []
        non_patho_mol_to_align = []
        patho_mol_to_align = []
        fasta_alignment = None
        non_pathogen_alignment = None
        pathogen_alignment = None
        spscore = None
        spscore_for_non_pathogen = None
        spscore_for_pathogen = None
        if glabrata_annotation.has_key('orthologs_in_candida_species'):
            for ortholog_in_candida in glabrata_annotation['orthologs_in_candida_species']:
                aligned_orthologs.append(ortholog_in_candida)
                ortholog_annotation = client[ortholog_in_candida.split('@')[1]]['annotations'].find_one({'_id': ortholog_in_candida.split('@')[0]})
                molecule_to_align = DNA(name=ortholog_annotation['locus_tag'], sequence=ortholog_annotation['translation'])# we get Candida orthologous sequence
                cgd_molecules_to_align.append(molecule_to_align)
                patho_mol_to_align.append(molecule_to_align)#all CGD Candida species are pathogenic
        if glabrata_annotation.has_key('sace_ortholog'):
            sace_molecule_to_align = DNA(name=glabrata_annotation['sace_ortholog'][0], sequence=glabrata_annotation['sace_ortholog'][1])# we get S. cerevisiae orthologous sequence
            cgd_molecules_to_align.append(sace_molecule_to_align)
            non_patho_mol_to_align.append(sace_molecule_to_align)

        ###  ALIGNMENT OF THE ORTHOLOGOUS PROTEIN SEQUENCES OF CANDIDA SPECIES & CEREVISIAE ###
        if cgd_molecules_to_align:
            aligned_orthologs.append("%s@Candida_glabrata_CBS_138"%glabrata_annotation['_id'])
            cagl_molecule_to_align = DNA(name=locus_tag, sequence=glabrata_annotation['translation'])# we get C. glabrata sequence
            cgd_molecules_to_align.append(cagl_molecule_to_align)
            patho_mol_to_align.append(cagl_molecule_to_align)
            muscle = Muscle(fasta_sequences=to_fasta(molecules=cgd_molecules_to_align))
            fasta_alignment = muscle.align()

            ### FOR EACH PROTEIN ALIGNMENT, SEACHING OF ORTHOLOGS IN THE NAKASEOMYCES SPECIES AND REALIGNMENT ###
            for species in nakaseo_dbs:
                hmmer = Hmmer(seqdb=nakaseo_dict[species][0]) # seqdb = BLAST protein database
                hmmer.hmmbuild(fasta_alignment, locus_tag)
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
                            if percent >= hmm_threshold:
                                for nakaseo_molecule in nakaseo_dict[species][1]:
                                    if nakaseo_molecule.name == data_hmmer['target']:
                                        if nakaseo_molecule.sequence.endswith('*'):
                                            nakaseo_molecule.sequence = nakaseo_molecule.sequence[:-1]
                                        nakaseo_molecules_to_align.append(nakaseo_molecule)
                                        if species in non_pathogens:
                                            non_patho_mol_to_align.append(nakaseo_molecule)
                                        else:
                                            patho_mol_to_align.append(nakaseo_molecule)
                                        aligned_orthologs.append("%s@%s"%(nakaseo_molecule.id, species))
                                        break
            
            ### REALIGNMENT WITH ALL ORTHOLOGOUS SEQUENCES ###
            if nakaseo_molecules_to_align:
                muscle = Muscle(fasta_sequences=to_fasta(molecules=nakaseo_molecules_to_align))
                if len(nakaseo_molecules_to_align) == 1:
                    fasta_alignment = muscle.realign(existing_msa=fasta_alignment)
                else:
                    input_msa = muscle.align()
                    fasta_alignment = muscle.realign(existing_msa=fasta_alignment, new_seqs=input_msa)

            ### SP SCORE CALCULATION FROM THE MSA ###
            muscle2 = Muscle(fasta_sequences=fasta_alignment)# MSA with or without Nakaseo
            spscore = muscle2.get_spscore()
            all_spscores.append(spscore)

            ### ALIGNMENT BETWEEN SEQUENCES OF NON-PATHOGENIC OR PATHOGENIC SPECIES ###
            if non_patho_mol_to_align:
                muscle = Muscle(fasta_sequences=to_fasta(molecules=non_patho_mol_to_align))
                non_pathogen_alignment = muscle.align()
                muscle2 = Muscle(fasta_sequences=non_pathogen_alignment)
                spscore_for_non_pathogen = muscle2.get_spscore()
                non_patho_spscores.append(spscore_for_non_pathogen)
            if patho_mol_to_align:
                muscle = Muscle(fasta_sequences=to_fasta(molecules=patho_mol_to_align))
                pathogen_alignment = muscle.align()
                muscle2 = Muscle(fasta_sequences=pathogen_alignment)
                spscore_for_pathogen = muscle2.get_spscore()
                patho_spscores.append(spscore_for_pathogen)

        ### UPDATING THE MONGO DATABASES ###
        if aligned_orthologs:
            alignment_id = str(ObjectId())
            for aligned_ortholog in aligned_orthologs:
                aligned_orthologs_copy = aligned_orthologs[:]
                aligned_orthologs_copy.remove(aligned_ortholog) # we remove the current ortholog in the list
                annotation = client[aligned_ortholog.split('@')[1]]['annotations'].find_one({'_id': aligned_ortholog.split('@')[0]})
                aln_ids = annotation.get('alignment')
                if not aln_ids:
                    client[aligned_ortholog.split('@')[1]]['annotations'].update({'_id': annotation['_id']},{'$set':{'alignment': [alignment_id], 'orthologs_in_candida_species': aligned_orthologs_copy}}, False)
                else:
                    ### A sequence that is orthologous to several genes of C. glabrata will be in several multiple alignments ### Example: CAAL locus "C4_02340W_A" ; CAGL loci "CAGL0M05665g" AND "CAGL0G06006g"
                    aln_ids.append(alignment_id)
                    candida_orthologs = annotation.get('orthologs_in_candida_species')
                    for ortho in aligned_orthologs_copy:
                        if ortho not in candida_orthologs:
                            candida_orthologs.append(ortho)
                    ### UPDATING THE CANDIDA AND NAKASEO MONGO DATABASES ###
                    client[aligned_ortholog.split('@')[1]]['annotations'].update({'_id': annotation['_id']},{'$set':{'alignment': aln_ids, 'orthologs_in_candida_species': candida_orthologs}}, False)

            ### INSERT A DOCUMENT IN THE COMPARATIVE GENOMICS MONGO DATABASE ###
            proteins_document = {
                                '_id': alignment_id,
                                'locus_tag': locus_tag,
                                'alignment': {'all_species': fasta_alignment},# if we need a msa at the format clustalw: convert_msa_format(fasta_alignment, "fasta", "clustal")
                                'spscore': {'all_species': spscore}
                                }

            if spscore_for_non_pathogen.__class__ is float:# we eliminate None but conserve 0.0 ; if spscore is not None, msa also
                proteins_document['alignment']['non_pathogenic_species'] = non_pathogen_alignment
                proteins_document['spscore']['non_pathogenic_species'] = spscore_for_non_pathogen
            if spscore_for_pathogen.__class__ is float:
                proteins_document['alignment']['pathogenic_species'] = pathogen_alignment
                proteins_document['spscore']['pathogenic_species'] = spscore_for_pathogen
            client['comparative_genomics']['proteins'].insert(proteins_document)

            purge_dir(directory="/tmp")
    purge_dir(directory="/tmp", option="blastdb")

    ### WE CALCULATE THE PERCENTILE OF SP SCORE FOR EACH ALIGNMENT ###
    for annotation in client['comparative_genomics']['proteins'].find({'spscore':{'$exists':True}}, no_cursor_timeout = True):
        percentile = {'all_species': scipy.stats.percentileofscore(all_spscores, annotation['spscore']['all_species'])}
        if annotation['spscore'].has_key('non_pathogenic_species'):
            percentile['non_pathogenic_species'] = scipy.stats.percentileofscore(all_spscores, annotation['spscore']['non_pathogenic_species'])
        if annotation['spscore'].has_key('pathogenic_species'):
            percentile['pathogenic_species'] = scipy.stats.percentileofscore(all_spscores, annotation['spscore']['pathogenic_species'])    
        client['comparative_genomics']['proteins'].update({'_id': annotation['_id']}, {'$set':{'percentile_of_spscore': percentile}}, False)

    client.close()

if __name__ == '__main__':
    hmm_threshold = 100

    if "-threshold" in sys.argv:
        hmm_threshold = int(sys.argv[sys.argv.index("-threshold")+1])

    hmm(hmm_threshold=hmm_threshold) 
