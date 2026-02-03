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
import time
import threading

try:
    from utils import shuffle_decoy_proteins, shuffle, shufflewithmut, read_fasta_file, digest, revswitch,writeseq, all_sublists, TRYPSIN,splitStringWithPeptide, get_new_peptide, get_new_pep_after_checkSimilar,concat_targe_decoy_protein,fasta_files_dedup
except:
    pass

try:
    import tqdm
except:
    print('Cannot import tqdm.py. Please install tqdm.')
    tqdm = False

def getDecoyProteinByRevert(args):
    '''get decoy proteins by revert
    '''
    
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
                if args.target_i2l:
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
    return upeps, dpeps


def shuffleForRevert(args, upeps, dpeps, peps_to_alt = None):
    
    # can only report total number in normal memory mode
    # find intersecting peptides
    if peps_to_alt is None:
        nonDecoys = upeps.intersection(dpeps)
    else:
        nonDecoys = peps_to_alt

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
            while ((aPep in upeps) or aPep == dPep) and i < args.maxit:

                # increment iteration counter
                i += 1

                # shuffle peptide
                aPep = shuffle(dPep, fix_C=args.cpos)

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
        return dAlternative
    else:
        os.rename(args.tout, args.dout)
        print("final decoy peptides:" + str(len(dpeps)))
        return {}

def get_target_peptides(args):
    '''digest proteins from input files. use when checkSimilar is enabled
    '''
    start_time = time.time()
    target_peptide_extra = set()
    for file_fasta in args.fasta:
        for header, seq in read_fasta_file(file_path=file_fasta):
            target_protein_seq = seq
            if not args.iso:
                target_protein_seq = target_protein_seq.replace('I', 'L')
            target_peptide_extra.update(TRYPSIN(target_protein_seq, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos))

    # add more upeps
    upeps_extra = target_peptide_extra
    print('target peptides after allowing missed clevage:',len(upeps_extra))
    print('checkSimilar is enabled!', args.checkSimilar)
    # if checkSimilar, N=D, Q=E,GG=N, replace N with D, replace Q with E, replace GG with N
    upeps_extra2 = set([get_new_pep_after_checkSimilar(i, args.checkSimilar) for i in upeps_extra])
    print('target peptides after', args.checkSimilar,len(upeps_extra2))
    del upeps_extra
    print("--- {:.2f} seconds digest input proteins allowing miss cleavage ---".format(time.time() - start_time))
    # open(args.tout + '.tempfile.upeps_extra2','w').write(str(upeps_extra2))
    args.upeps_extra2 = upeps_extra2

def get_decoy_peptides(args):
    '''digest proteins from final decoy file. use when checkSimilar is enabled
    '''
    start_time = time.time()
    decoy_peptide_extra = set()
    for header, seq in read_fasta_file(file_path=args.dout):
        decoy_protein_seq = seq
        if not args.iso:
            decoy_protein_seq = decoy_protein_seq.replace('I', 'L')
        decoy_peptide_extra.update(TRYPSIN(decoy_protein_seq, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos))

    # add more upeps
    dpeps_extra = decoy_peptide_extra
    print('decoy peptides after allowing missed clevage:',len(dpeps_extra))
    print('checkSimilar is enabled!',args.checkSimilar)
    # if checkSimilar, N=D, Q=E,GG=N, replace N with D, replace Q with E, replace GG with N
    dpeps_extra2 = set([get_new_pep_after_checkSimilar(i, args.checkSimilar) for i in dpeps_extra])
    print('decoy peptides after ', args.checkSimilar, len(dpeps_extra2))
    del dpeps_extra
    print("--- {:.2f} seconds digest decoy proteins allowing miss cleavage ---".format(time.time() - start_time))
    # open(args.tout + '.tempfile.upeps_extra2','w').write(str(upeps_extra2))
    args.dpeps_extra2 = dpeps_extra2

def get_decoy_proteins(args, fout):
    '''use when checkSimilar is enabled. 
    get decoy proteins with peptides in upeps
    '''
    start_time = time.time()
    # get those do not need further change and save. 
    ls_decoy_proteins = []#store those need change in ls_decoy_proteins
    upeps_extra2 = args.upeps_extra2
    n_no_change = 0
    # add more upeps
    decoy_peptide_extra = set()
    for header, seq in read_fasta_file(args.tout):
        if args.iso:
            seq_for_digest = seq
        else:
            seq_for_digest = seq.replace('I', 'L')
        ls_decoy_pep = TRYPSIN(seq_for_digest, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)
        decoy_peptide_extra.update(ls_decoy_pep)
        ls_decoy_tochange = [i for i in ls_decoy_pep if get_new_pep_after_checkSimilar(i, args.checkSimilar) in upeps_extra2]
        if len(ls_decoy_tochange) == 0:
            n_no_change += 1
            fout.write(header+'\n'+seq+'\n')
        else:
            ls_decoy_proteins.append([header, seq, set(ls_decoy_tochange)])
    
    dpeps_extra = decoy_peptide_extra
    print('decoy peptides after allowing missed clevage:',len(dpeps_extra))
    dpeps_extra = set([ get_new_pep_after_checkSimilar(i, args.checkSimilar) for i in dpeps_extra])
    print('decoy peptides after', args.checkSimilar, len(dpeps_extra))
    print(n_no_change, 'decoy proteins do not contain similar peptide in target proteins')
    print('number of decoy peptides to alter',sum(len(i[2]) for i in ls_decoy_proteins), 'of which', len(set(i for j in ls_decoy_proteins for i in j[2])),'are unique')
    print('number of proteins to alter', len(ls_decoy_proteins))
    del dpeps_extra
    print("--- {:.2f} seconds get decoy proteins with peptides in target proteins ---".format(time.time() - start_time))

    # open(args.tout + '.tempfile.ls_decoy_proteins','w').write(str(ls_decoy_proteins))
    return ls_decoy_proteins


def get_amino_acid_frequency(args):
    '''get amino acid frequency in file
    '''
    start_time = time.time()
    amino_acids = Counter()
    for header, seq in read_fasta_file(args.tout):
        amino_acids += Counter(seq)
    ## exclude non normal amino acids
    normal_AA = 'ADEFGHLMSTVWYRKNQPCI'
    amino_acids = {k:v for k,v in amino_acids.items() if k in normal_AA}
    amino_acids = {k:v/sum(amino_acids.values()) for k,v in amino_acids.items()}
    print('amino acid frequencies in input proteins are:', {k:str(round(v,3)) for k,v in amino_acids.items()})
    if not args.iso:
        print('Note: I were changed to L')
    args.amino_acids = amino_acids
    print("--- {:.2f} seconds amino acid frequency ---".format(time.time() - start_time))


def update_one_protein(header, seq, ls_decoy_tochange, dAlternative2, ls_peptide_altered, upeps_extra2, args):
    n_decoy_tochange = len(ls_decoy_tochange)
    random.shuffle(ls_decoy_tochange)
    iii = 0
    while(len(ls_decoy_tochange) != 0 and iii <= 2 * n_decoy_tochange):
        iii += 1
        for peptide in ls_decoy_tochange:
            solved_by_dAlternative2 = False
            
            # first try to alter with existing peptides in dAlternative2
            if peptide in dAlternative2:
                for new_pep in dAlternative2[peptide]:
                    l = splitStringWithPeptide(proseq = seq, peptide = peptide, anti_cleavage_sites = args.noc, cleavage_sites = args.csites)
                    if len(ls_decoy_tochange) == 1 and len(l) == 1:
                        print('warning!', header, seq, ls_decoy_tochange, l)
                    l_new = [new_pep if i == peptide else i for i in l]
                    proseq_changed_by_dAlternative = ''.join(l_new)
                    ls_decoy_pep = TRYPSIN(proseq_changed_by_dAlternative, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)
                    ls_decoy_pep = list(set([i for i in ls_decoy_pep if get_new_pep_after_checkSimilar(i, args.checkSimilar) in upeps_extra2]))
                    if len(ls_decoy_pep) < len(ls_decoy_tochange):
                        # keep the change
                        solved_by_dAlternative2 = True
                        seq = proseq_changed_by_dAlternative
                        ls_decoy_tochange = ls_decoy_pep
                        break
            # if not in dAlternative2, or cannot get good result by using the one in dAlternative2, create new peptide
            if (peptide not in dAlternative2) or (solved_by_dAlternative2 == False):
                l = splitStringWithPeptide(proseq = seq, peptide = peptide, anti_cleavage_sites = args.noc, cleavage_sites = args.csites)
                if len(ls_decoy_tochange) == 1 and len(l) == 1:
                    print('warning!', header, seq, ls_decoy_tochange, l)
                for i in range(4 * args.maxit):
                    pep_mut_start_step = i // args.maxit # start from different step 
                    new_pep, new_pep_comment = get_new_peptide(peptide, upeps_extra2, start_step = pep_mut_start_step,indel_ratio = 1, amino_acids = args.amino_acids, fix_C = args.cpos, maxit=100,checkSimilar=args.checkSimilar)
                    l_new = [new_pep if i == peptide else i for i in l]
                    proseq_changed_by_dAlternative = ''.join(l_new)
                    ls_decoy_pep = TRYPSIN(proseq_changed_by_dAlternative, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)
                    ls_decoy_pep = list(set([i for i in ls_decoy_pep if get_new_pep_after_checkSimilar(i, args.checkSimilar) in upeps_extra2]))
                    if len(ls_decoy_pep) < len(ls_decoy_tochange):
                        # keep the change
                        solved_by_dAlternative2 = True
                        seq = proseq_changed_by_dAlternative
                        ls_decoy_tochange = ls_decoy_pep
                        if peptide not in dAlternative2:
                            dAlternative2[peptide] = []
                        if new_pep not in dAlternative2[peptide]:
                            dAlternative2[peptide].append(new_pep)
                            ls_peptide_altered.append((peptide, new_pep, new_pep_comment))
                            # print((peptide, new_pep, new_pep_comment))
                        break
            if solved_by_dAlternative2:
                break
    
    return header, seq, ls_decoy_tochange


def checkSimilarForProteins(args):
    '''check similar proteins and change
    '''
    os.rename(args.dout, args.tout)
    dAlternative = args.dAlternative

    # digest input files allowing miss cleavage
    thread_upeps = threading.Thread(target = get_target_peptides, args = (args,))
    thread_upeps.start()
    
    # get amino acid frequency
    thread_aa = threading.Thread(target = get_amino_acid_frequency, args = (args,))
    thread_aa.start()

    thread_upeps.join()
    thread_aa.join()
    upeps_extra2 = args.upeps_extra2
    amino_acids = args.amino_acids

    # save decoy proteins with no overlapped peptides in upeps_extra2 in args.tout. for the rest, save ls_decoy_proteins in file args.tout + '.temp.pickle'
    fout = open(args.dout,'w')
    ls_decoy_proteins = get_decoy_proteins(args, fout)# count of decoy proteins to alter
    n_decoy_proteins = len(ls_decoy_proteins)
    args.n_decoy_proteins = n_decoy_proteins


    # try to solve shuffle peptide only
    ## update dAlternative2
    print('try to solve by shuffle peptide')
    dAlternative_similar = {k:v for k,v in dAlternative.items() if v != '' and get_new_pep_after_checkSimilar(v, args.checkSimilar) not in upeps_extra2}

    # shuffle create a new dict to store dAlternative_similar. use list as values for multiple choices of alternative peptides
    dAlternative2 = {k:[v] for k,v in dAlternative_similar.items()}

    # alt proteins
    ls_peptide_altered = []
    ls_unsolved_proteins = []
    n_improved = 0
    n_solved = 0
    n_noimprove = 0
    if tqdm:
        iter_ls_decoy_proteins = tqdm.tqdm(range(len(ls_decoy_proteins)))
    else:
        iter_ls_decoy_proteins = range(len(ls_decoy_proteins))
    for n in iter_ls_decoy_proteins:
        header, seq, ls_decoy_tochange = ls_decoy_proteins[n]
        ls_decoy_tochange = list(set(ls_decoy_tochange))
        n_decoy_tochange = len(ls_decoy_tochange)
        header, seq2, ls_decoy_tochange2 = update_one_protein(header, seq, ls_decoy_tochange, dAlternative2, ls_peptide_altered, upeps_extra2, args)
        
        ls_decoy_pep = TRYPSIN(seq2, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)
        ls_decoy_pep = list(set([i for i in ls_decoy_pep if get_new_pep_after_checkSimilar(i, args.checkSimilar) in upeps_extra2]))
        if ls_decoy_pep != ls_decoy_tochange2 or len(ls_decoy_pep) != 0:
            print('warning, the {n}th protein cannot be solved by shuffling peptides. sequence is {seq}, altered seq is {seq2}. oringal decoy peptides: {ls_decoy_pep}. after shuffling, the decoy peptides: {ls_decoy_tochange2}'.format(n=n,ls_decoy_tochange2=ls_decoy_tochange2,ls_decoy_pep=ls_decoy_pep, seq=seq, seq2=seq2))
            # break
        
        # save the protein sequence
        fout.write(header+'\n'+seq2+'\n')
        if len(ls_decoy_tochange2) == 0:
            n_solved += 1
        elif len(ls_decoy_tochange2) < n_decoy_tochange:
            n_improved += 1
            ls_unsolved_proteins.append([header, seq2, ls_decoy_tochange2])
        else:
            n_noimprove += 1
            ls_unsolved_proteins.append([header, seq2, ls_decoy_tochange2])


    print('number of unsolved decoy peptides:', len(set([j for i in ls_unsolved_proteins for j in i[2]])))
    print('number of proteins with unsolved decoy peptides:', n_improved + n_noimprove)
    fout.close()

    # delete temporary file
    os.remove(args.tout)

    get_decoy_peptides(args)
    print('double check decoy proteins. number of overlapped peptides considering', args.checkSimilar, len(args.dpeps_extra2 & upeps_extra2))


def main():
    args.tout = args.dout + '.tempfile'
    if not args.target_file:
        args.target_file = args.dout + '.temptargetfile'

    checkSimilar = args.checkSimilar
    if checkSimilar:
        print('run orginal DecoyPYrat pipeline first. checkSimilar is enabled; I/L are treated as the same unless --no_isobaric is set; shuffle is enabled')
        args.noshuf = False # have to shuffle
    
    start_time = time.time()

    if args.dedup:
        print('dedup is enabled. First, remove duplicated sequences in the input file.')
        fasta_files_dedup(file_input = args.fasta, file_output = args.dout + '.temp_targetfile_dedup', ILdifferent = (not args.target_i2l))
        args.fasta = [args.dout + '.temp_targetfile_dedup']
    

    if not args.all_shuffle_mimic:
        upeps, dpeps = getDecoyProteinByRevert(args)
        dAlternative = shuffleForRevert(args, upeps, dpeps)
        del upeps, dpeps
    else:
        print('all_shuffle_mimic is enabled. Get decoy proteins by revert first. Then all decoy peptides are randomized')
        args.noshuf = False # have to shuffle
        upeps, dpeps = getDecoyProteinByRevert(args)
        dAlternative = shuffleForRevert(args, upeps, dpeps,peps_to_alt=dpeps)
        del upeps, dpeps
    print("--- {:.2f} seconds run DecoyPYrat pipeline. Note: pseudo-reverse method was used in default.  ---".format(time.time() - start_time))

    if checkSimilar:
        args.dAlternative = dAlternative
        checkSimilarForProteins(args)
    
    if args.concat != '':
        print('concat is set to {}. This setting is best for multistage search'.format(args.concat))
        os.rename(args.dout, args.tout)
        concat_targe_decoy_protein(file_target = args.target_file, file_decoy = args.tout, file_output = args.dout, concat=args.concat, keepname = args.names)
    
    # remove temporary files
    if args.target_file == args.dout + '.temptargetfile':
        os.remove(args.target_file)
    if os.path.exists(args.tout):
        os.remove(args.tout)
    if os.path.exists(args.dout + '.temp_targetfile_dedup'):
        os.remove(args.dout + '.temp_targetfile_dedup')
    
    print("--- {:.2f} seconds run the whole pipeline ---".format(time.time() - start_time))



if __name__ == "__main__":
    # Read command line arguments and create help documentation using argparse
    parser = argparse.ArgumentParser(
            description='''Create decoy protein sequences. Each protein is reversed and the cleavage sites switched with preceding amino acid. 
            Peptides are checked for existence in target sequences if found the tool will attempt to shuffle them.
            ''')
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
                        help='Turn OFF shuffling of decoy peptides that are in the target database. Default=False')
    parser.add_argument('--all_shuffle_mimic', '-R', dest='all_shuffle_mimic', default=False, action='store_true',
                        help='use random seq when generating decoy proteins. Similar method like mimic. Default=False')
    parser.add_argument('--do_not_switch', '-s', dest='noswitch', default=False, action='store_true',
                        help='''If set, revert the whole sequence as decoy. Default=False. Use the pseudo-reverse method used by Comet.''')
    parser.add_argument('--decoy_prefix', '-d', dest='dprefix', default='DECOY',
                        help='Set accesion prefix for decoy proteins in output. Default=DECOY')
    parser.add_argument('--output_fasta', '-o', dest='dout', default='decoy.fa',
                        help='Set file to write decoy proteins to. Default=decoy.fa')
    parser.add_argument('--no_isobaric', '-i', dest='iso', default=False,
                        action='store_true', help='Do not make decoy peptides isobaric. Default=False, I will be changed to L in decoy sequences')
    parser.add_argument('--target_I2L', dest='target_i2l', default=False,
                        action='store_true', help='Convert I to L in target output and dedup. Default=False')
    parser.add_argument('--threads', '-N', dest='threads', default=1, type=int,
                        help='number of threads to use. default 1. Note: currently only one thread is allowed')
    parser.add_argument('--keep_names', '-k', dest='names', default=False,
                        action='store_true', help='''Keep sequence names in the decoy output. Default=False. if decoy_prefix = "DECOY", name will be like DECOY_1, DECOY_2. The orginal name will be included in the header line, separated by "\t\t". If set, name will be like DECOY_{oringal_name}''')
    parser.add_argument('--target', '-T', dest='target_file', default="",
                        help='Combine and store the target file. Target sequences are unchanged by default; use --target_I2L to convert I to L. Default="", do not save file')
    parser.add_argument('--checkSimilar', '-S', dest='checkSimilar', default=False, 
                        help='''If set, I/L are treated as the same unless --no_isobaric is set; target output I to L is controlled by --target_I2L.
                        To enable this option, it is recommended to set checkSimilar="GG=N,N=D,Q=E".
                        Allow missed cleavage, and max_peptide_length will be used when determining shared peptides in target and decoy sequences. In default setting, the digested peptides do not overlap with each other. 
                        "Peptides are checked for existence in target sequences and if found the tool will attempt to shuffle them iterativly until they are unique". 
                        Additionally, consider those amino acids were equal: N=D, Q=E, GG=N. 
                        If cannot solve after max_iterations of shuffling, introduce AA mutations, deletions or insertions. AA not in cleavage_sites and anti_cleavage_sites. The number_of_changes (sum of mutations, deletions and insertions) <= 1 for each peptide. 
                        2 * max_iterations times.
                        do_not_shuffle option will be ignored.
                        Default=False''')
    parser.add_argument('--concat','-C', dest='concat', default='', help='''default = "", do not append target sequences to decoy sequences. If concat is set, the target sequence is appended to decoy sequences, with concat as the separator. use '*' for comet or other engines which treat "*" as separator. the same as cleavage_sites ("R" in most cases) for other search engines ''')
    parser.add_argument('--dedup','-U', dest='dedup', default=False, action='store_true', help='''If set, try to remove duplicated sequences in the target database. The first sequence will be kept. default = False.''')
    args = parser.parse_args()
    # args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz  /data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar GG=N,N=D,Q=E --target target.fa '''.split())
    # args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz  /data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar GG=N,N=D,Q=E --target target.fa --threads 10'''.split())
    # args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz  /data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar GG=N,N=D,Q=E --target target.fa -R'''.split())
    main()
