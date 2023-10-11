
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
import random
import time
import gc

def read_fasta_file(file_path):
    """
    Reads a fasta file and header and sequence.
    header include '>', but not '\n'
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
        yield header, sequence.strip('*')
    
    file.close()

# Tryptic Digest - Can be modified to take 'sites' as argument and digest based on that
def digest(protein, sites='KR', pos='c', no='P', min_len=0):
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
    l = list(filter(lambda x: len(x) >= min_len, (protein.split(','))))
    return [i for i in l if i]

def get_new_pep_after_checkSimilar(seq, checkSimilar='GG=N,N=D,Q=E'):
    '''
    do replacement of Amino Acids as described in checkSimilar. In default setting, replace GG to N, then N to D, then Q to E.
    '''
    seqnew = seq
    l = checkSimilar.split(',')
    for e in l:
        old, new = e.split('=')
        seqnew = seqnew.replace(old, new)
    return seqnew

def split_into_n_parts(lst, n):
    """Split a list into n equal parts"""
    part_size = len(lst) // n
    remainder = len(lst) % n
    
    parts = []
    for i in range(n):
        if i < remainder:
            parts.append(lst[i*part_size + i:(i+1)*part_size + i + 1])
        else:
            parts.append(lst[i*part_size + remainder:(i+1)*part_size + remainder])
    return parts

def split_into_n_parts_equal_step(lst, n):
    """Split a list into n equal parts"""
    parts = [[] for i in range(n)]
    for i in range(len(lst)):
        parts[i % n].append(i)
    return parts


def revswitch(protein, noswitch=False, csites='KR', noc='P', cpos='c'):
    
    """Return a reversed protein sequence with cleavage residues switched with preceding residue
    protein is the protein sequence
    if noswitch == True, revert the sequence. Else, use the pseudo-reversed method like Comet
    csites, cleavage sites
    noc, anti cleavage sites
    cpos, c or n, cleavage at the n or c terminal
    """
    # reverse protein sequence with a reverse splice convert to list
    if noswitch:
        return protein[::-1]

    if noswitch == False:
        l = digest(protein=protein, sites=csites, pos=cpos, no=noc, min_len=0)
        for n in range(len(l)):
            e = l[n]
            if cpos == 'c':
                if e[-1] in csites:
                    e = e[:-1][::-1] + e[-1]
                else:
                    e = e[::-1]
            elif cpos == 'n':
                if e[0] in csites:
                    e = e[0] + e[1:][::-1]
                else:
                    e = e[::-1]
            else:
                print('error, csites must be n or c')
            l[n] = e
    
        # return reversed with/without switched proteins as string
        return ''.join(l)


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
        s = ''
    random.shuffle(l)
    # return new peptide
    return ''.join(l) + s


def shufflewithmut(peptide, indel_ratio = 0.1, amino_acids = None, fix_C = True):
    """shuffle peptide without moving c-terminal amino acid cleavage site, with one mutation, insertion or deletion
    if indel_ratio > 1, it will be treated as the count of mutation, insertions or deletions
    """
    if fix_C:
        # extract terminal aa
        s = peptide[-1]
        # convert peptide to list (remove K/R) and shuffle the list
        l = list(peptide[:-1])
        rand_pos = random.randint(0, len(l) - 1)
    else:
        # extract terminal aa
        s = ''
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

def shufflewithmutMultiple(peptide, indel_ratio = 0.1, amino_acids = None, fix_C = True, count_mut_sites = 1):
    '''same as above, but introducing count_mut_sites mutations
    '''
    for i in range(count_mut_sites):
        new_pep = shufflewithmut(peptide, indel_ratio, amino_acids, fix_C)
    return new_pep


def get_new_peptide(peptide, upeps_extra2, start_step = 0, indel_ratio = 1, amino_acids = None, fix_C = True, maxit=100,checkSimilar='GG=N,N=D,Q=E'):
    '''shuffle peptide until get a peptide not in upeps_extra2. 4 steps.
    first, shuffle with function shuffle, for max maxit rounds
    secondly, shuffle with shufflewithmut, for max maxit rounds, indel_ratio set to 0
    thirdly, shuffle with shufflewithmut, for max maxit rounds, indel_ratio set to 1
    fourth, shuffle with shufflewithmutMultiple, for max maxit rounds, indel_ratio set to 0.5, count_mut_sites set to 2
    start_step = [0,1,2]. start from 0, 1, or 2 steps
    '''
    if start_step == 0:
        for i in range(maxit):
            new_pep = shuffle(peptide)
            if get_new_pep_after_checkSimilar(new_pep, checkSimilar) not in upeps_extra2:
                return new_pep,'shuffle'
    
    if start_step <= 1:
        for i in range(maxit):
            new_pep = shufflewithmut(peptide, indel_ratio=0, amino_acids=amino_acids,fix_C=fix_C)
            if get_new_pep_after_checkSimilar(new_pep, checkSimilar) not in upeps_extra2:
                return new_pep,'mut'
    
    if start_step <= 2:
        for i in range(maxit):
            new_pep = shufflewithmut(peptide, indel_ratio=1, amino_acids=amino_acids,fix_C=fix_C)
            if get_new_pep_after_checkSimilar(new_pep, checkSimilar) not in upeps_extra2:
                return new_pep,'indel'
    
    if start_step <= 3:
        for i in range(maxit):
            new_pep = shufflewithmutMultiple(peptide, indel_ratio=0.5, amino_acids=amino_acids,fix_C=fix_C, count_mut_sites=2)
            if get_new_pep_after_checkSimilar(new_pep, checkSimilar) not in upeps_extra2:
                return new_pep,'mutIndel'
    
    return peptide, 'noAlternative'


def writeseq(args, seq, upeps, dpeps, outfa, pid, dcount):
    # make sequence isobaric (check args for switch off)
    if args.iso == False:
        seq = seq.replace('I', 'L')

    # digest sequence add peptides to set
    upeps.update(digest(seq, args.csites, args.cpos, args.noc, args.minlen))

    # reverse and switch protein sequence
    decoyseq = revswitch(seq, args.noswitch, args.csites, args.noc, args.cpos)

    # update decoy peptide set
    dpeps.update(digest(decoyseq, args.csites,
                    args.cpos, args.noc, args.minlen))

    # write decoy protein accession and sequence to file
    if args.names:
        outfa.write('>{}_{}\n'.format(args.dprefix, pid.strip(">")))
    else:
        outfa.write('>' + args.dprefix + '_' + str(dcount) + '\t\t' + pid.strip(">") +'\n')
    outfa.write(decoyseq + '\n')

def all_sublists(lst):
    '''return a sublist from lst
    '''
    for i in range(len(lst) +1):
        for sublist in itertools.combinations(lst,i):
            yield(sublist)




def TRYPSIN(protein, sites, pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40):
    '''
    Digests a protein sequence using trypsin and returns the resulting peptides within a specified length range.

    Args:
    protein (str): The amino acid sequence of the protein to be digested.
    sites (str): A string of amino acid cleavage sites for trypsin.
    pos (str, optional): A string indicating whether to perform n-terminal or c-terminal cleavage (default is 'c').
    no (str, optional): Amino acids following the cleavage site that would prevent cleavage (default is 'P').
    miss_cleavage (int, optional): The maximum number of missed cleavage sites allowed (default is 2).
    peplen_min (int, optional): The minimum length of peptides to be returned (default is 6).
    peplen_max (int, optional): The maximum length of peptides to be returned (default is 40).

    Returns:
    A list of peptides resulting from trypsin digestion of the protein sequence within the specified length range.
   
    '''
    peptides_cut_all = digest(protein = protein, sites=sites, pos=pos, no=no, min_len=0)
    peptides = []
    for i in range(miss_cleavage + 1):
        for j in range(len(peptides_cut_all) - i):
            peptide = ''.join(peptides_cut_all[j:j + i + 1])
            peplen = len(peptide)
            if peplen >= peplen_min and peplen <= peplen_max:
                peptides.append(peptide)

    return peptides


def splitStringWithPeptide(proseq, peptide, anti_cleavage_sites='P', cleavage_sites='KR',pos='c'):
    '''split proseq to parts, separated by peptide
    anti_cleavage_sites is the AA after the peptide, like "P"
    cleavage_sites is the AA before the peptide, like "KR"
    pos is n or c. cut at n or c terminal
    '''
    if not peptide:
        return [proseq]
    if peptide not in proseq:
        return [proseq]
    
    # count number of miss cleavage in peptide
    pep_from_pep = digest(protein = peptide, sites=cleavage_sites, pos=pos, no=anti_cleavage_sites, min_len=0)
    pep_from_pep = [i for i in pep_from_pep if i] # remove empty string
    n_miss_cleavage_site = len(pep_from_pep) - 1
    if n_miss_cleavage_site > 10:
        print(proseq, peptide,'double check cleavage')
    
    # digest protein
    peptides = digest(protein = proseq, sites=cleavage_sites, pos=pos, no=anti_cleavage_sites, min_len=0)
    peptides = [i for i in peptides if i]# remove empty string
    positions = [[0, len(peptides[0])]]
    for p in peptides[1:]:
        positions.append([positions[-1][1], positions[-1][1] + len(p)])
    # get peptides with n_miss_cleavage_site
    peptides_miss_cleavage = []
    positions_miss_cleavage = []
    positions_peptide = []
    for i in range(len(peptides) - n_miss_cleavage_site):
        pep_miss = ''.join(peptides[i:i + n_miss_cleavage_site + 1])
        pep_miss_pos = [positions[i][0], positions[i+ n_miss_cleavage_site][1]]
        peptides_miss_cleavage.append(pep_miss)
        positions_miss_cleavage.append(pep_miss_pos)
        if pep_miss == peptide:
            positions_peptide.append(pep_miss_pos)
    
    # remove overlapped positions
    while True:
        for i in range(len(positions_peptide) - 1):
            if positions_peptide[i][1] > positions_peptide[i+1][0]:
                positions_peptide.pop(i + 1)
                break
        else:
            break
    
    break_points = [0] + [i for j in positions_peptide for i in j] + [len(proseq)]
    positions_peptide = [[break_points[i],break_points[i+1]] for i in range(len(break_points) -1)]
    peptides_seg = [proseq[i:j] for i,j in positions_peptide]
    return [i for i in peptides_seg if i]


def get_new_protein_with_pep_mut_multiple(ls_decoy_proteins,new_decoy_peptides, dAlternative2, alter_protein_better,upeps_extra2, args):
    results = []
    for header, seq, ls_decoy_tochange in ls_decoy_proteins:
        ls_pep_to_change = set([p for p in ls_decoy_tochange if p in dAlternative2])
        if len(ls_pep_to_change) == 0:
            results.append([[header, seq, ls_decoy_tochange],[], False])
            continue
        ls_decoy_tochange = sorted(ls_decoy_tochange, key=lambda x:new_decoy_peptides[x] if x in new_decoy_peptides else float('inf'), reverse=True) # change the most abundant peptides first. if peptide not in new_decoy_peptides, it should be put in front as the peptide was identified before
        peptide_changed = []
        for p in ls_pep_to_change:
            l = splitStringWithPeptide(proseq = seq, peptide = p, anti_cleavage_sites = args.noc, cleavage_sites = args.csites)
            if len(l) > 1 and p in dAlternative2:
                peptide_changed.append(p)
            if len(ls_decoy_tochange) == 1 and len(l) == 1:
                print('warning!', header, seq, ls_decoy_tochange, l)
            l = [dAlternative2[i] if i in dAlternative2 else i for i in l]
            proseq_changed_by_dAlternative = ''.join(l)
            
        ls_decoy_pep = TRYPSIN(proseq_changed_by_dAlternative, miss_cleavage=args.miss_cleavage, peplen_min=args.minlen, peplen_max=args.maxlen, sites=args.csites, no=args.noc, pos=args.cpos)
        ls_decoy_tochange2 = set([i for i in ls_decoy_pep if get_new_pep_after_checkSimilar(i, args.checkSimilar) in upeps_extra2])
        if alter_protein_better:
            good = len(ls_decoy_tochange2) < len(ls_decoy_tochange)
        else:
            good = len(ls_decoy_tochange2) <= len(ls_decoy_tochange)
        if good:
            results.append([[header, proseq_changed_by_dAlternative, ls_decoy_tochange2], peptide_changed, True])
        else:
            results.append([[header, seq, ls_decoy_tochange],[], False])
    
    return results

def get_new_protein_with_pep_mut_multiple_disk(file_ls_decoy_proteins, alter_protein_better, args, ls_decoy_proteins_ranges = None):
    ls_decoy_proteins = eval(open(file_ls_decoy_proteins).read())
    if ls_decoy_proteins_ranges is not None:
        ls_decoy_proteins = [ls_decoy_proteins[i] for i in ls_decoy_proteins_ranges]
    
    new_decoy_peptides = eval(open(args.tout + '.tempfile.new_decoy_peptides').read())
    dAlternative2 = eval(open(args.tout + '.tempfile.dAlternative2').read())
    upeps_extra2 = eval(open(args.tout + '.tempfile.upeps_extra2').read())
    return get_new_protein_with_pep_mut_multiple(ls_decoy_proteins,new_decoy_peptides, dAlternative2, alter_protein_better,upeps_extra2, args)

def update_dAlternative2(args, shuffle_method = 'shuffle',indel_ratio=0,amino_acids=None,fix_C=True):
    '''update dAlternative2 by shuffling
    '''
    file_ls_decoy_proteins = args.tout + '.tempfile.ls_decoy_proteins'
    ls_decoy_proteins = eval(open(file_ls_decoy_proteins).read())
    file_dAlternative2 = args.tout + '.tempfile.dAlternative_similar'
    dAlternative2 = eval(open(file_dAlternative2).read())
    file_upeps_extra2 = args.tout + '.tempfile.upeps_extra2'
    upeps_extra2 = eval(open(file_upeps_extra2).read())

    new_decoy_peptides = [i for j in ls_decoy_proteins for i in j[2]]
    n_pep_count_all = len(new_decoy_peptides) # all counts of decoy peptides to mute
    new_decoy_peptides = Counter(new_decoy_peptides)
    n_pep_count_unique = len(new_decoy_peptides)# unique new decoy peptides to mutate
    new_decoy_peptides = {k:v for k,v in new_decoy_peptides.items() if k not in dAlternative2}
    n_pep_existing = n_pep_count_unique - len(new_decoy_peptides)
    n_pep_shuffle_new = 0
    n_pep_shufflemut_new = 0
    n_pep_cannot_solve = 0
    for p in new_decoy_peptides:
        if len(set(p)) == 1 and shuffle_method =='shuffle':
            n_pep_cannot_solve += 1
        else:
            for i in range(args.maxit):
                if shuffle_method =='shuffle':
                    new_pep = shuffle(p)
                elif shuffle_method == 'shufflewithmut':
                    new_pep = shufflewithmut(p, indel_ratio=indel_ratio, amino_acids=amino_acids,fix_C=fix_C)
                if get_new_pep_after_checkSimilar(new_pep, args.checkSimilar) not in upeps_extra2:
                    dAlternative2[p] = new_pep
                    if shuffle_method == 'shuffle':
                        n_pep_shuffle_new += 1
                    elif shuffle_method =='shufflewithmut':
                        n_pep_shufflemut_new += 1
                    break
            else:
                n_pep_cannot_solve += 1


    txt = f'total peptides to alt: {n_pep_count_all}, unique peptides count is: {n_pep_count_unique}, of them, '
    if n_pep_shuffle_new != 0:
        txt += f'peptide changed by shuffle: {n_pep_shuffle_new}; '
    if n_pep_shufflemut_new != 0:
        txt += f'peptide changed by shuffle with one mutation: {n_pep_shufflemut_new}; '
    if n_pep_cannot_solve != 0:
        txt += f'peptide with no alternative choices: {n_pep_cannot_solve}; '
    if n_pep_existing != 0:
        txt += f'peptide existed previously: {n_pep_existing}'
    print(txt)
    file_dAlternative2 = args.tout + '.tempfile.dAlternative2'
    file_new_decoy_peptides = args.tout + '.tempfile.new_decoy_peptides'
    open(file_dAlternative2,'w').write(str(dAlternative2))
    open(file_new_decoy_peptides, 'w').write(str(new_decoy_peptides))
    # return dAlternative2, new_decoy_peptides


def split_ls_decoy_proteins_to_n_parts(file_ls_decoy_proteins, nparts):
    '''file_ls_decoy_proteins stores ls_decoy_proteins
    '''
    ls_decoy_proteins = eval(open(file_ls_decoy_proteins).read())
    ls_decoy_proteins = sorted(ls_decoy_proteins, key=lambda x:len(set(x[2])), reverse=True)
    ls_ranges = split_into_n_parts_equal_step(range(len(ls_decoy_proteins)), nparts)
    ls_file_ls_decoy_proteins = []
    for n in range(nparts):
        f = file_ls_decoy_proteins + '.part.' + str(n)
        ls_file_ls_decoy_proteins.append(f)
        open(f,'w').write(str([ls_decoy_proteins[i] for i in ls_ranges[n]]))
    
    return ls_file_ls_decoy_proteins


def shuffle_decoy_proteins(fout, args, shuffle_method = 'shuffle', **kwargs):
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
    if 'alter_protein_better' in kwargs:
        alter_protein_better = kwargs['alter_protein_better']
    else:
        alter_protein_better = True
    if 'threads' in kwargs:
        threads = kwargs['threads']
    
    n_decoy_proteins = args.n_decoy_proteins
    ## get new dAlternative2 by shuffling peptides
    update_dAlternative2(shuffle_method = shuffle_method, args=args,indel_ratio=indel_ratio,amino_acids=amino_acids,fix_C=fix_C)
    ## change protein sequences and check
    file_ls_decoy_proteins = args.tout + '.tempfile.ls_decoy_proteins'
    if args.threads == 1 or n_decoy_proteins < 10000:
        results = get_new_protein_with_pep_mut_multiple_disk(file_ls_decoy_proteins, alter_protein_better, args, ls_decoy_proteins_ranges = None)
    else:
        nparts = max(args.threads, n_decoy_proteins // 10000) #split to args.threads parts or len(ls_decoy_proteins) // 10000 parts, whichever is larger
        ls_file_ls_decoy_proteins = split_ls_decoy_proteins_to_n_parts(file_ls_decoy_proteins, nparts)
        ls_params = [[i, alter_protein_better, args] for i in ls_file_ls_decoy_proteins]
        # random.shuffle(ls_decoy_proteins)# avoid long running time for some parts
        pool = Pool(args.threads)
        # pool = multiprocessing.pool.ThreadPool(args.threads)
        results = pool.starmap(get_new_protein_with_pep_mut_multiple_disk, ls_params, chunksize=1)
        pool.close()
        pool.join()
        for f in ls_file_ls_decoy_proteins:
            os.remove(f)
        os.remove(file_ls_decoy_proteins)
        results = [i for j in results for i in j]
    ls_decoy_proteins = [i[0] for i in results]
    ## save proteins
    for n in range(len(ls_decoy_proteins)):
        header, seq, ls_decoy_tochange = ls_decoy_proteins[n]
        if len(ls_decoy_tochange) == 0:
            fout.write(header+'\n'+seq+'\n')
    ls_decoy_proteins = [i for i in ls_decoy_proteins if len(i[2]) != 0]
    n_pep_shuffle_accepted = sum(len(i[1]) for i in results)
    n_pr_changed = sum(i[2] for i in results)
    print('number of total proteins: {}, number proteins changed: {}, number of proteins left: {}\nnumber of shuffled peptides acepted: {}, of them, {} were unique '.format(n_decoy_proteins,n_pr_changed,len(ls_decoy_proteins),n_pep_shuffle_accepted, len(set([i for j in results for i in j[1]]))))
    open(args.tout + '.tempfile.ls_decoy_proteins','w').write(str(ls_decoy_proteins))
    return len(ls_decoy_proteins)


def concat_targe_decoy_protein(file_target, file_decoy, file_output, concat = '*', keepname = False):
    '''combine sequences in file_target and file_decoy, save in file_output. sequences joined by concat.
    if keepname == True, sequences
    '''
    dc_target = {}
    for header, seq in read_fasta_file(file_path=file_target):
        seq = seq.upper().strip('*')
        dc_target[header] = seq
    
    fout_decoy = open(file_output,'w')
    for header, seq in read_fasta_file(file_path=file_decoy):
        seq = seq.upper().strip('*')
        if keepname:
            seq_new = seq + concat + dc_target['>' + header.split('_', maxsplit=1)[-1]]
        else:
            seq_new = seq + concat + dc_target['>' + header.split('\t\t')[-1]]
        fout_decoy.write('{}\n{}\n'.format(header, seq_new))
    fout_decoy.close()


def fasta_files_dedup(file_input, file_output, ILdifferent = False):
    '''remove duplicated sequences in file_input, save in file_output.
    if ILdifferent == True, do not change I to L.
    else, change I to L.
    '''
    if isinstance(file_input, str):
        file_input = [file_input]

    # open file to save the result
    outfa = open(file_output, 'w')

    n_count = 0
    dc_seqcount = {}
    # Open FASTA file using first cmd line argument
    for file_fasta in file_input:
        for header, seq in read_fasta_file(file_path=file_fasta):
            n_count += 1
            seq = seq.upper().strip('*')

            if ILdifferent == False:
                seq = seq.replace('I', 'L')
            
            if seq not in dc_seqcount:
                dc_seqcount[seq] = 1
                outfa.write('{}\n{}\n'.format(header, seq))
            else:
                dc_seqcount[seq] += 1


    outfa.close()
    print('number of input sequences: {}. number of unique sequences: {}'.format(n_count, len(dc_seqcount)))
    print('{sequence counts: number of sequences}:',dict(Counter(dc_seqcount.values())))