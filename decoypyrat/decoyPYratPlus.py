# DecoyPYrat - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses
#
# MIT License
#
# Copyright (c) 2016 James Christopher Wright - Wellcome Trust Sanger Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# used to get cmd line arguments
import argparse
# used to shuffle peptides
import random
# used to rename/delete tmp file
import os
import gzip
import itertools
from collections import Counter
import numpy as np
import re
from multiprocessing import Pool
import multiprocessing
import pickle

from utils import shuffle_decoy_proteins, shuffle, shufflewithmut, read_fasta_file, digest, revswitch,writeseq, all_sublists, TRYPSIN,splitStringWithPeptide


def main():
    checkSimilar = args.checkSimilar
    if checkSimilar:
        print('checkSimilar is enabled, I is always replaced with L.')
        args.iso =False
    

    # Create empty sets to add all target and decoy peptides
    upeps = set()
    dpeps = set()

    # Counter for number of decoy sequences
    dcount = 1

    # empty protein sequence
    seq = ''

    # open temporary decoy FASTA file
    outfa = open(args.tout, 'w')

    
    if args.target_file:
        outfa_target = open(args.target_file, 'w')

    # Open FASTA file using first cmd line argument
    for file_fasta in args.fasta:
        for header, seq in read_fasta_file(file_path=file_fasta):
            seq = seq.upper().strip('*')
            writeseq(args, seq, upeps, dpeps, outfa, header, dcount)
            if args.target_file:
                if args.iso == False:
                    outfa_target.write('{}\n{}\n'.format(header, seq.replace('I', 'L')))
                else:
                    outfa_target.write('{}\n{}\n'.format(header, seq))
            dcount += 1
    
    outfa.close()

    if args.target_file:
        outfa_target.close()

    # Summarise the numbers of target and decoy peptides and their intersection
    nonDecoys = set()
    print("proteins:" + str(dcount - 1))
    print("target peptides:" + str(len(upeps)))
    print("decoy peptides:" + str(len(dpeps)))

    # run original command to generate decoy file
    if checkSimilar:
        print('first, run orginal DecoyPYrat pipeline')
    # Reloop decoy file in reduced memory mode to store only intersecting decoys

    # can only report total number in normal memory mode
    # find intersecting peptides
    nonDecoys = upeps.intersection(dpeps)

    print("#intersection:" + str(len(nonDecoys)))
    # if there are decoy peptides that are in the target peptide set
    if len(nonDecoys) > 0 and args.noshuf == False:

        # create empty dictionary with bad decoys as keys
        dAlternative = dict.fromkeys(nonDecoys, '')
        noAlternative = list()

        # loop bad decoys / dictionary keys
        for dPep in dAlternative:
            i = 0
            aPep = dPep

            # shuffle until aPep is not in target set (maximum of 10 iterations)
            while aPep in upeps and i < args.maxit:

                # increment iteration counter
                i += 1

                # shuffle peptide
                aPep = shuffle(dPep)

            # update dictionary with alternative shuffled peptide
            dAlternative[dPep] = aPep

            # warn if peptide has no suitable alternative, add to removal list
            if i == args.maxit:
                noAlternative.append(dPep)

        print(str(len(noAlternative)) + ' have no alternative peptide')
        # remove peptides with no alternative
        for p in noAlternative:
            del dAlternative[p]

        # open second decoy file
        with open(args.dout, "wt") as fout:
            # open original decoy file
            with open(args.tout, "rt") as fin:
                # loop each line of original decoy fasta
                for line in fin:
                    # if line is not accession replace peptides in dictionary with alternatives
                    if line[0] != '>':
                        # digest decoy sequence
                        for p in digest(line.rstrip(), args.csites, args.cpos, args.noc, args.minlen):
                            # store decoy peptide for final count
                            dpeps.add(p)

                            # if decoy peptide is in dictionary replace with alternative
                            if p in dAlternative:
                                line = line.replace(p, dAlternative[p])

                    fout.write(line)
            fin.close()
        fout.close()

        # delete temporary file
        os.remove(args.tout)
    else:
        os.rename(args.tout, args.dout)

        print("final decoy peptides:" + str(len(dpeps)))

    if checkSimilar:
        os.rename(args.dout, args.tout)

        target_peptide_extra = []
        for file_fasta in args.fasta:
            for header, seq in read_fasta_file(file_path=file_fasta):
                target_protein_seq = seq
                target_protein_seq = target_protein_seq.replace('I', 'L')
                target_peptide_extra += TRYPSIN(target_protein_seq, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)

        # add more upeps
        upeps_extra = set(target_peptide_extra)
        print('target peptides after allowing missed clevage:',len(upeps_extra))
        print('checkSimilar is enabled! N=D, Q=E, GG=N')
        # if checkSimilar, N=D, Q=E,GG=N, replace N with D, replace Q with E, replace GG with N
        upeps_extra2 = set([i.replace('N','D').replace('Q','E').replace('GG', 'N') for i in upeps_extra])
        print('target peptides after N=D, Q=E, GG=N:',len(upeps_extra2))
        del upeps_extra, upeps, dpeps

        
        # get those do not need further change and save. 
        ls_decoy_proteins = []#store those need change in ls_decoy_proteins
        fout = open(args.dout,'w')
        n_no_change = 0
        # add more upeps
        decoy_peptide_extra = []
        for header, seq in read_fasta_file(args.tout):
            ls_decoy_pep = TRYPSIN(seq.replace('I', 'L'), miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)
            decoy_peptide_extra += ls_decoy_pep
            ls_decoy_tochange = [i for i in ls_decoy_pep if i.replace('N','D').replace('Q','E').replace('GG', 'N') in upeps_extra2]
            if len(ls_decoy_tochange) == 0:
                n_no_change += 1
                fout.write(header+'\n'+seq+'\n')
            else:
                ls_decoy_proteins.append([header, seq, set(ls_decoy_tochange)])
        
        dpeps_extra = set(decoy_peptide_extra)
        print('decoy peptides after allowing missed clevage:',len(dpeps_extra))
        dpeps_extra = set([ i.replace('N','D').replace('Q','E').replace('GG', 'N') for i in dpeps_extra])
        print('decoy peptides after N=D, Q=E, GG=N:', len(dpeps_extra))
        del target_peptide_extra, decoy_peptide_extra
        print(n_no_change, 'decoy proteins do not contain similar peptide in target proteins')
        print('number of decoy peptides to alter:', len(set(i for j in ls_decoy_proteins for i in j[2])))
        print('number of proteins to alter', len(ls_decoy_proteins))

        # get amino acid frequency
        amino_acids = Counter()
        for header, seq in read_fasta_file(args.tout):
            amino_acids += Counter(seq)
        ## exclude non normal amino acids
        normal_AA = 'ADEFGHLMSTVWYRKNQPCI'
        amino_acids = {k:v for k,v in amino_acids.items() if k in normal_AA}
        amino_acids = {k:v/sum(amino_acids.values()) for k,v in amino_acids.items()}
        print('amino acid frequencies in input proteins are:', {k:str(round(v,3)) for k,v in amino_acids.items()})
        print('Note: I were changed to L')

        # try to solve shuffle peptide only
        ## update dAlternative2
        print('try to solve by shuffle peptide only.')
        dAlternative_similar = {k:v for k,v in dAlternative.items() if v != '' and v.replace('N','D').replace('Q','E').replace('GG', 'N') not in upeps_extra2}

        for nround in range(10):
            # create dAlternative2
            print(nround,'round of shuffling peptide')
            dAlternative2 = dAlternative_similar.copy()
            ls_decoy_proteins = shuffle_decoy_proteins(ls_decoy_proteins, dAlternative2, args=args, upeps_extra2=upeps_extra2,fout=fout, shuffle_method = 'shuffle')
            if len(ls_decoy_proteins) == 0:
                break

        ## shuffle with mutation allowed
        indel_rate = 0.1
        n = 2
        if args.all_shuffle_mimic:
            n = 5
        for nround in range(args.maxit * n):
            print(nround,'round of shuffling peptide, allow one mutation per peptide')
            # create dAlternative2
            dAlternative2 = dAlternative_similar.copy()
            if nround > 50 or len(ls_decoy_proteins) < 50:
                indel_rate += 0.002
                indel_rate = min(indel_rate, 0.8)
                print(f'increase mutation indel_ratio from 0.1 to {indel_rate:.3f}')
            ls_decoy_proteins = shuffle_decoy_proteins(ls_decoy_proteins, dAlternative2, args=args,amino_acids=amino_acids, upeps_extra2=upeps_extra2,fout=fout, shuffle_method = 'shufflewithmut', indel_ratio = indel_rate)
            if len(ls_decoy_proteins) == 0:
                break

        

        print('number of unsolved decoy peptides:', len(set([j for i in ls_decoy_proteins for j in i[2]])))
        print('number of proteins with unsolved decoy peptides:', len(ls_decoy_proteins))
        for header, seq, ls_decoy_tochange in ls_decoy_proteins:
            fout.write(header+'\n'+seq+'\n')
        fout.close()

        # delete temporary file
        os.remove(args.tout)



if __name__ == "__main__":
    # Read command line arguments and create help documentation using argparse
    parser = argparse.ArgumentParser(
            description='''Create decoy protein sequences. Each protein is reversed and the cleavage sites switched with preceding amino acid. 
            Peptides are checked for existence in target sequences if found the tool will attempt to shuffle them.
            James.Wright@sanger.ac.uk 2015''')
    parser.add_argument('fasta', metavar='*.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz',
                        nargs='+',
                        help='FASTA file of target proteins sequences for which to create decoys')
    parser.add_argument('--cleavage_sites', '-c', dest='csites', default='KR',
                        help='A list of amino acids at which to cleave during digestion. Default = KR')
    parser.add_argument('--anti_cleavage_sites', '-a', dest='noc', default='P',
                        help='A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = "P"')
    parser.add_argument('--cleavage_position', '-p', dest='cpos', default='c', choices=[
                        'c', 'n'], help='Set cleavage to be c or n terminal of specified cleavage sites. Default = c')
    parser.add_argument('--min_peptide_length', '-l', dest='minlen', default=6, type=int,
                        help='Set minimum length of peptides to compare between target and decoy. Default = 6')
    parser.add_argument('--max_peptide_length', '-M', dest='maxlen', default=40, type=int,
                        help='Set max length of peptides to compare between target and decoy.  Default = 40')
    parser.add_argument('--max_iterations', '-n', dest='maxit', default=100, type=int,
                        help='Set maximum number of times to shuffle a peptide to make it non-target before failing. Default=100')
    parser.add_argument('--miss_cleavage', '-L', dest='miss_cleavage', default=2, type=int,
                        help='miss_cleavage when digesting protein. only work when checkSimilar is enabled. Default=2')
    parser.add_argument('--do_not_shuffle', '-x', dest='noshuf', default=False, action='store_true',
                        help='Turn OFF shuffling of decoy peptides that are in the target database. Default=false')
    parser.add_argument('--all_shuffle_mimic', '-R', dest='all_shuffle_mimic', default=False, action='store_true',
                        help='use random seq when generating decoy proteins. Similar method like mimic. Default=false')
    parser.add_argument('--do_not_switch', '-s', dest='noswitch', default=False, action='store_true',
                        help='Turn OFF switching of cleavage site with preceding amino acid. Default=false')
    parser.add_argument('--decoy_prefix', '-d', dest='dprefix', default='XXX',
                        help='Set accesion prefix for decoy proteins in output. Default=XXX')
    parser.add_argument('--output_fasta', '-o', dest='dout', default='decoy.fa',
                        help='Set file to write decoy proteins to. Default=decoy.fa')
    parser.add_argument('--temp_file', '-t', dest='tout', default='tmp.fa',
                        help='Set temporary file to write decoys prior to shuffling. Default=tmp.fa')
    parser.add_argument('--no_isobaric', '-i', dest='iso', default=False,
                        action='store_true', help='Do not make decoy peptides isobaric. Default=false, I will be changed to L in decoy sequences')
    parser.add_argument('--threads', '-N', dest='threads', default=1, type=int,
                        help='number of threads to use. default 1')
    parser.add_argument('--keep_names', '-k', dest='names', default=False,
                        action='store_true', help='Keep sequence names in the decoy output. Default=false')
    parser.add_argument('--target', '-T', dest='target_file', default="",
                        help='Combine and store the target file. I will be changed to L default. If no_isobaric, I will not be changed. Default="", do not save file')
    parser.add_argument('--checkSimilar', '-S', dest='checkSimilar', default=False, action='store_true',
                        help='''If set, ItoL is enabled automatically; output_fasta will include target sequences by changing I to L; allow overlapped digestion, and max_peptide_length will be used. In default setting, the digested peptides do not overlap with each other. 
                        "Peptides are checked for existence in target sequences and if found the tool will attempt to shuffle them iterativly until they are unique". 
                        Additionally, consider those amino acids were equal: N=D, Q=E, GG=N. 
                        If cannot solve after max_iterations of shuffling, introduce AA mutations, deletions or insertions. AA not in cleavage_sites and anti_cleavage_sites. The number_of_changes (sum of mutations, deletions and insertions) <= 1 for each peptide. 
                        2 * max_iterations times.
                        do_not_shuffle option will be ignored.
                        Default=false''')
    args = parser.parse_args()
    # args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz  /data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar --target target.fa '''.split())
    # args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz  /data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar --target target.fa --threads 10'''.split())
    # args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz  /data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar --target target.fa -R'''.split())
    main()
