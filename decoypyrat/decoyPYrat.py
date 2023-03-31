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
parser.add_argument('--anti_cleavage_sites', '-a', dest='noc', default='',
                    help='A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = none')
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
parser.add_argument('--memory_save', '-m', dest='mem', default=False, action='store_true',
                    help='Slower but uses less memory (does not store decoy peptide list). Default=false')
parser.add_argument('--keep_names', '-k', dest='names', default=False,
                    action='store_true', help='Keep sequence names in the decoy output. Default=false')
parser.add_argument('--target', '-T', dest='target_file', default="",
                    help='Combine and store the target file. I will be changed to L default. If no_isobaric, I will not be changed. Default="", do not save file')
parser.add_argument('--checkSimilar', '-S', dest='checkSimilar', default=False, action='store_true',
                    help='''If set, ItoL is enabled automatically; output_fasta will include target sequences by changing I to L; allow overlapped digestion, and max_peptide_length will be used. In default setting, the digested peptides do not overlap with each other. 
                    "Peptides are checked for existence in target sequences and if found the tool will attempt to shuffle them iterativly until they are unique". 
                    Additionally, consider those amino acids were equal: N=D, Q=E, GG=N. 
                    If cannot solve after max_iterations of shuffling, introduce AA mutations, deletions or insertions. AA not in cleavage_sites and anti_cleavage_sites. The number_of_changes (sum of mutations, deletions and insertions) <= 1 for each peptide. 
                     3 * max_iterations times.
                     do_not_shuffle option will be ignored.
                     Default=false''')

args = parser.parse_args('''/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090.fasta.gz
/data/pub/genome/mouse/uniprot/20230301/UP000000589_10090_additional.fasta.gz /data/p/maxquant/MaxQuant_2.1.4.0/bin/conf/contaminants.fasta --checkSimilar --target target.fa
'''.split())


def read_fasta_file(file_path):
    """
    Reads a fasta file and header and sequence.
    """
    if file_path.endswith('.gz'):
        file = gzip.open(file_path,'rt')
    else:
        file = open(file_path)
    header = ''
    sequence = ''
    for line in file:
        if line.startswith('>'):
            if sequence:
                yield header, sequence
                sequence = ''
            header = line.strip()
        else:
            sequence += line.strip()
    if sequence:
        yield header, sequence
    
    file.close()

# Tryptic Digest - Can be modified to take 'sites' as argument and digest based on that
def digest(protein, sites, pos, no, min):
    """Return a list of cleaved peptides with minimum length in protein sequence.
            protein = sequence
            sites = string of amino acid cleavage sites
            pos = n or c for n-terminal or c-terminal cleavage
            no = amino acids following site that would prevent cleavage ie proline
            min = minimum length of peptides returned"""

    # for each possible cleavage site insert a comma with before or after depending on pos
    for s in sites:
        r = s + ','
        if pos == 'n':
            r = ',' + s
        protein = protein.replace(s, r)

    # for each possible cleavage and all none cleavage remove comma
    for s in sites:
        for n in no:
            a = s + ',' + n
            if pos == 'n':
                a = ',' + s + n
            r = s + n
            protein = protein.replace(a, r)

    # filter peptides into list by minimum size
    return list(filter(lambda x: len(x) >= min, (protein.split(','))))


def revswitch(protein, noswitch, sites):
    """Return a reversed protein sequence with cleavage residues switched with preceding residue"""
    # reverse protein sequence with a reverse splice convert to list
    revseq = list(protein[::-1])

    if noswitch == False:

        # loop sequence list
        for i, c in enumerate(revseq):
            # if value is cleavage site switch with previous amino acid
            for s in sites:
                if c == s:
                    aa = revseq[i-1]
                    revseq[i-1] = revseq[i]
                    revseq[i] = aa

    # return reversed with/without switched proteins as string
    return ''.join(revseq)


def shuffle(peptide, fix_C = True):
    """shuffle peptide without moving c-terminal amino acid cleavage site if fix_C == True.
    if fix_C == False, shuffle the whole peptide
    """
    if fix_C:
        # extract terminal aa
        s = peptide[-1]
        # convert peptide to list (remove K/R) and shuffle the list
        l = list(peptide[:-1])
    else:
        l = list(peptide)
    random.shuffle(l)
    # return new peptide
    return ''.join(l) + s


def shufflewithmut(peptide, indel_ratio = 0.1, amino_acids = None, fix_C = True):
    """shuffle peptide without moving c-terminal amino acid cleavage site, with one mutation, insertion or deletion"""
    if fix_C:
        # extract terminal aa
        s = peptide[-1]
        # convert peptide to list (remove K/R) and shuffle the list
        l = list(peptide[:-1])
        rand_pos = random.randint(0, len(l) - 1)
    else:
        # extract terminal aa
        s = peptide
        # convert peptide to list (remove K/R) and shuffle the list
        l = list(peptide)
        rand_pos = random.randint(0, len(l))
    if amino_acids is None:
        amino_acids = list('ADEFGHLMSTVWYRKNQPCI')# 20 AA
        AA_freq = [1/len(amino_acids) for _ in amino_acids]
    else:
        AA_freq = list(amino_acids.values())
        amino_acids = list(amino_acids.keys())
    
    new_AA = np.random.choice(amino_acids, p=AA_freq)
    if random.random() < indel_ratio:
        if random.random() < 0.5:
            l[rand_pos] += new_AA #introduce insertion in 10% chance
        else:
            l[rand_pos] = '' # introduce deletion
    else:
        l[rand_pos] = new_AA
    random.shuffle(l)
    # return new peptide
    return ''.join(l) + s

# def shuffleProteinWithmut(protein, anti_cleavage_sites, cleavage_sites, mut = 0, indel_ratio = 0.1, amino_acids = None):
#     """shuffle protein without change position of cleavage_sites and anti_cleavage_sites
#     if mut is the number of mutations. introduce with mutation, insertion or deletion randomly"""
#     l = list(protein)
#     if amino_acids is None:
#         amino_acids = list('ADEFGHLMSTVWYRKNQPCI')# 20 AA
    
#     for i in range(mut):
#         rand_pos = random.randint(0, len(l) - 1)
#         new_AA = random.choice(amino_acids)
#         if random.random() < indel_ratio:
#             if random.random() < 0.5:
#                 l[rand_pos] += new_AA #introduce insertion in 10% chance
#             else:
#                 l[rand_pos] = '' # introduce deletion
#         else:
#             l[rand_pos] = new_AA

#     protein = ''.join(l)
    
#     l = list(protein)
#     l_fix = []
#     l_nonfix = []
#     for i in range(len(protein)):
#         if protein[i] in anti_cleavage_sites or protein[i] in cleavage_sites:
#             l_fix.append(i)
#         else:
#             l_nonfix.append(i)
    
#     random.shuffle(l_nonfix)
#     l_new = []
#     for i in range(len(protein)):
#         if i in l_fix:
#             l_new.append(protein[i])
#         else:
#             l_new.append(protein[l_nonfix[-1]])
#             l_nonfix.pop()

#     # return new peptide
#     return ''.join(l_new)

def writeseq(args, seq, upeps, dpeps, outfa, pid, dcount):
    # make sequence isobaric (check args for switch off)
    if args.iso == False:
        seq = seq.replace('I', 'L')

    # digest sequence add peptides to set
    upeps.update(digest(seq, args.csites, args.cpos, args.noc, args.minlen))

    # reverse and switch protein sequence
    decoyseq = revswitch(seq, args.noswitch, args.csites)

    # do not store decoy peptide set in reduced memory mode
    if args.mem == False:
        # update decoy peptide set
        dpeps.update(digest(decoyseq, args.csites,
                     args.cpos, args.noc, args.minlen))

    # write decoy protein accession and sequence to file
    if args.names:
        outfa.write('>{}_{}\n'.format(args.dprefix, pid))
    else:
        outfa.write('>' + args.dprefix + '_' + str(dcount) + '\n')
    outfa.write(decoyseq + '\n')

def all_sublists(lst):
    '''return a sublist from lst
    '''
    for i in range(len(lst) +1):
        for sublist in itertools.combinations(lst,i):
            yield(sublist)


def shuffle_decoy_proteins(ls_decoy_proteins, dAlternative2, fout, args, upeps_extra2, shuffle_method = 'shuffle', **kwargs):
    '''
    amino_acids = None, fix_C = True
    '''
    if 'indel_ratio' in kwargs:
        indel_ratio = kwargs['indel_ratio']
    else:
        indel_ratio = 0.1
    if 'amino_acids' in kwargs:
        amino_acids = kwargs['amino_acids']
    else:
        amino_acids = None
    if 'fix_C' in kwargs:
        fix_C = kwargs['fix_C']
    else:
        fix_C = True
    
    n_pep_shuffle_new = 0
    n_pep_shufflemut_new = 0
    n_pep_cannot_solve = 0
    # update dAlternative2
    for header, seq, ls_decoy_tochange in ls_decoy_proteins:
        for p in ls_decoy_tochange:
            if p not in dAlternative2:
                for i in range(args.maxit):
                    if shuffle_method =='shuffle':
                        new_pep = shuffle(p)
                    elif shuffle_method == 'shufflewithmut':
                        new_pep = shufflewithmut(p, indel_ratio=indel_ratio, amino_acids=amino_acids,fix_C=fix_C)
                    if new_pep.replace('N','D').replace('Q','E').replace('GG', 'N') not in upeps_extra2:
                        dAlternative2[p] = new_pep
                        if shuffle_method == 'shuffle':
                            n_pep_shuffle_new += 1
                        elif shuffle_method =='shufflewithmut':
                            n_pep_shufflemut_new += 1
                        break
                else:
                    n_pep_cannot_solve += 1
    n = n_pep_shuffle_new + n_pep_shufflemut_new + n_pep_cannot_solve
    print(f'total peptides to alt: {n}; peptide changed by shuffle: {n_pep_shuffle_new}; peptide changed by shuffle with one mutation: {n_pep_shufflemut_new}; peptide with no alternative choices: {n_pep_cannot_solve}')
    
    ## change protein sequences and check
    for n in range(len(ls_decoy_proteins)):
        header, seq, ls_decoy_tochange = ls_decoy_proteins[n]
        ls_decoy_tochange = list(ls_decoy_tochange)
        random.shuffle(ls_decoy_tochange)
        for p in ls_decoy_tochange:
            l = splitStringWithPeptide(proseq = seq, peptide = p, anti_cleavage_sites = args.noc, cleavage_sites = args.csites)
            l = [dAlternative2[i] if i in dAlternative2 else i for i in l]
            proseq_changed_by_dAlternative = ''.join(l)
        ls_decoy_pep = TRYPSIN(proseq_changed_by_dAlternative, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, csites=args.csites, noc=args.noc)
        ls_decoy_tochange2 = set([i for i in ls_decoy_pep if i.replace('N','D').replace('Q','E').replace('GG', 'N') in upeps_extra2])
        if len(ls_decoy_tochange2) < len(ls_decoy_tochange):
            ls_decoy_proteins[n] = [header, proseq_changed_by_dAlternative, ls_decoy_tochange2]

    ## save proteins
    for n in range(len(ls_decoy_proteins)):
        header, seq, ls_decoy_tochange = ls_decoy_proteins[n]
        if len(ls_decoy_tochange) == 0:
            fout.write(header+'\n'+seq+'\n')
    ls_decoy_proteins = [i for i in ls_decoy_proteins if len(i[2]) != 0]
    print('number of proteins left', len(ls_decoy_proteins))
    return ls_decoy_proteins


def TRYPSIN(proseq, miss_cleavage=2, peplen_min=6, peplen_max=40, csites='KR', noc='P'):
    '''
    code modified from
    https://github.com/yafeng/trypsin
    only return peptides >= peplen_min and <= peplen_max
    csites is cleavage site, default is "KR"
    noc is non-cleavage site, default is "P"
    '''
    peptides = []
    cut_sites = [0]
    for i in range(0, len(proseq)-1):
        if proseq[i] in csites and proseq[i+1] not in noc:
            cut_sites.append(i+1)

    if cut_sites[-1] != len(proseq):
        cut_sites.append(len(proseq))

    number_of_cut_sites = len(cut_sites)

    for cut_start, cut_end in itertools.combinations(range(number_of_cut_sites), 2):
        if cut_end - cut_start <= miss_cleavage + 1:
            peplen = cut_sites[cut_end] - cut_sites[cut_start]
            if peplen >= peplen_min and peplen <= peplen_max:
                peptide = proseq[cut_sites[cut_start]:cut_sites[cut_end]]
                peptides.append(peptide)

    return peptides


def splitStringWithPeptide(proseq, peptide, anti_cleavage_sites, cleavage_sites):
    '''split proseq to parts, separated by peptide
    anti_cleavage_sites is the AA after the peptide, like "P"
    cleavage_sites is the AA before the peptide, like "KR"
    '''
    if not peptide:
        return [proseq]

    i = 0
    positions = []
    while i <= len(proseq) - len(peptide):
        pep_start = proseq.find(peptide, i)
        if pep_start < 0:
            break
        elif pep_start == 0:
            pep_end = pep_start + len(peptide)
            if pep_end == len(proseq):
                positions.append([pep_start, pep_end])
            elif proseq[pep_end] not in anti_cleavage_sites:
                positions.append([pep_start, pep_end])
                i = pep_end
            else:
                i += 1
        else:
            pep_end = pep_start + len(peptide)
            if proseq[pep_start -1] not in cleavage_sites:
                i += 1
            elif pep_end == len(proseq):
                positions.append([pep_start, pep_end])
            elif proseq[pep_end] not in anti_cleavage_sites:
                positions.append([pep_start, pep_end])
                i = pep_end
            else:
                i += 1

    break_points = [0] + [i for j in positions for i in j] + [len(proseq)]
    positions = [[break_points[i],break_points[i+1]] for i in range(len(break_points) -1)]
    return [proseq[i:j] for i,j in positions]

def main():
    # Create empty sets to add all target and decoy peptides
    upeps = set()
    dpeps = set()

    # Counter for number of decoy sequences
    dcount = 1

    # empty protein sequence
    seq = ''

    # open temporary decoy FASTA file
    outfa = open(args.tout, 'w')

    checkSimilar = args.checkSimilar
    if args.target_file:
        outfa_target = open(args.target_file, 'w')
    
    if checkSimilar:
        target_peptide_extra = []

    # Open FASTA file using first cmd line argument
    for file_fasta in args.fasta:
        for header, seq in read_fasta_file(file_path=file_fasta):
            seq = seq.upper()
            writeseq(args, seq, upeps, dpeps, outfa, header, dcount)
            if args.target_file:
                if args.iso == False:
                    outfa_target.write('{}\n{}\n'.format(header, seq.replace('I', 'L')))
                else:
                    outfa_target.write('{}\n{}\n'.format(header, seq))
            dcount += 1
            if checkSimilar:
                target_protein_seq = seq
                if args.iso == False:
                    target_protein_seq = target_protein_seq.replace('I', 'L')
                target_peptide_extra += TRYPSIN(target_protein_seq, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, csites=args.csites, noc=args.noc)
    
    outfa.close()

    if args.target_file:
        outfa_target.close()

    # Summarise the numbers of target and decoy peptides and their intersection
    nonDecoys = set()
    print("proteins:" + str(dcount - 1))
    print("target peptides:" + str(len(upeps)))
    print("decoy peptides:" + str(len(dpeps)))

    # add more upeps
    if checkSimilar:
        upeps_extra = upeps | set(target_peptide_extra)
        print('target peptides after allowing missed clevage:',len(upeps_extra))
        # add more upeps
        decoy_peptide_extra = []
        for header, seq in read_fasta_file(args.tout):
            decoy_peptide_extra += TRYPSIN(seq, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, csites=args.csites, noc=args.noc)

        dpeps_extra = dpeps | set(decoy_peptide_extra)
        print('decoy peptides after allowing missed clevage:',len(dpeps_extra))
        del target_peptide_extra, decoy_peptide_extra
    
    # run original command to generate decoy file
    if checkSimilar:
        print('first, run orginal DecoyPYrat pipeline')
    # Reloop decoy file in reduced memory mode to store only intersecting decoys
    if args.mem:
        # open temp decoys
        with open(args.tout, "rt") as fin:
            for line in fin:
                # if line is not accession
                if line[0] != '>':
                    # digest protein
                    peps = digest(line.rstrip(), args.csites, args.cpos, args.noc, args.minlen)
                    for p in p:
                        # check if in target peptides if true then add to nonDecoys
                        if p in upeps:
                            nonDecoys.add(p)
        fin.close()
        print("decoy peptides: !Memory Saving Made!")
    else:
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
        # if checkSimilar, N=D, Q=E,GG=N, replace N with D, replace Q with E, replace GG with N
        print('checkSimilar is enabled! N=D, Q=E, GG=N')
        upeps_extra2 = set([i.replace('N','D').replace('Q','E').replace('GG', 'N') for i in upeps_extra])
        print('target peptides after N=D, Q=E, GG=N:',len(upeps_extra2))

        # get those do not need further change and save. 
        ls_decoy_proteins = []#store those need change in ls_decoy_proteins
        fout = open(args.dout,'w')
        n_no_change = 0
        for header, seq in read_fasta_file(args.tout):
            ls_decoy_pep = TRYPSIN(seq, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, csites=args.csites, noc=args.noc)
            ls_decoy_tochange = [i for i in ls_decoy_pep if i.replace('N','D').replace('Q','E').replace('GG', 'N') in upeps_extra2]
            if len(ls_decoy_tochange) == 0:
                n_no_change += 1
                fout.write(header+'\n'+seq+'\n')
            else:
                ls_decoy_proteins.append([header, seq, set(ls_decoy_tochange)])
        
        print(n_no_change, 'decoy proteins do not contain similar peptide in target proteins')

        # get amino acid frequency
        amino_acids = Counter()
        for header, seq in read_fasta_file(args.tout):
            amino_acids += Counter(seq)
        ## exclude non normal amino acids
        normal_AA = 'ADEFGHLMSTVWYRKNQPCI'
        amino_acids = {k:v for k,v in amino_acids.items() if k in normal_AA}
        amino_acids = {k:v/sum(amino_acids.values()) for k,v in amino_acids.items()}

        # try to solve shuffle peptide only
        ## update dAlternative2
        print('try to solve by shuffle peptide only.')
        for nround in range(10):
            # create dAlternative2
            print(nround,'round of shuffling peptide')
            dAlternative2 = {k:v for k,v in dAlternative.items() if v != '' and v.replace('N','D').replace('Q','E').replace('GG', 'N') not in upeps_extra2}
            ls_decoy_proteins = shuffle_decoy_proteins(ls_decoy_proteins, dAlternative2, args=args, upeps_extra2=upeps_extra2,fout=fout, shuffle_method = 'shuffle')
            if len(ls_decoy_proteins) == 0:
                break

        ## shuffle with mutation allowed
        for nround in range(args.maxit):
            print(nround,'round of shuffling peptide, allow one mutation per peptide')
            # create dAlternative2
            dAlternative2 = {k:v for k,v in dAlternative.items() if v != '' and v.replace('N','D').replace('Q','E').replace('GG', 'N') not in upeps_extra2}
            ls_decoy_proteins = shuffle_decoy_proteins(ls_decoy_proteins, dAlternative2, args=args,amino_acids=amino_acids, upeps_extra2=upeps_extra2,fout=fout, shuffle_method = 'shufflewithmut', indel_ratio = 0.1)
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
    main()
