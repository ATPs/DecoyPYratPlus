import argparse
import gzip
import os
import sys
import time
from functools import lru_cache

_FAST_DIGEST_ENABLED = False
_FAST_DIGEST_AVAILABLE = False
_FAST_DIGEST_MODULE = None


def _load_fast_digest():
    try:
        import fast_digest as fast_mod
        return fast_mod
    except Exception:
        pass
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(base_dir)
        if parent_dir and parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)
        try:
            import fast_digest as fast_mod
            return fast_mod
        except Exception:
            from DecoyPYratPlus import fast_digest as fast_mod
            return fast_mod
    except Exception:
        return None


def set_fast_digest(enabled):
    """Enable or disable fast digest/trypsin (Cython). Falls back to Python if unavailable."""
    global _FAST_DIGEST_ENABLED, _FAST_DIGEST_AVAILABLE, _FAST_DIGEST_MODULE
    _FAST_DIGEST_ENABLED = bool(enabled)
    if not _FAST_DIGEST_ENABLED:
        return False
    if _FAST_DIGEST_AVAILABLE:
        return True
    fast_mod = _load_fast_digest()
    if fast_mod is None:
        try:
            import pyximport
            pyximport.install(language_level=3)
            fast_mod = _load_fast_digest()
        except Exception:
            fast_mod = None
    if fast_mod is None or not hasattr(fast_mod, "digest") or not hasattr(fast_mod, "trypsin"):
        _FAST_DIGEST_ENABLED = False
        _FAST_DIGEST_AVAILABLE = False
        _FAST_DIGEST_MODULE = None
        print("fast_digest unavailable; falling back to pure Python.")
        return False
    _FAST_DIGEST_AVAILABLE = True
    _FAST_DIGEST_MODULE = fast_mod
    return True


def read_fasta_file(file_path):
    """
    Reads a fasta file and yields header and sequence.
    header includes '>', but not '\n'
    """
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt') as file_handle:
        header = ''
        sequence_parts = []
        for line in file_handle:
            if line.startswith('>'):
                if sequence_parts:
                    yield header, ''.join(sequence_parts)
                    sequence_parts = []
                header = line.strip()
            else:
                sequence_parts.append(line.strip())
        if sequence_parts:
            yield header, ''.join(sequence_parts).strip('*')


def _digest_python(protein, sites='KR', pos='c', no='P', min_len=0):
    """Return a list of cleaved peptides with minimum length in protein sequence."""
    for s in sites:
        replacement = s + ','
        if pos == 'n':
            replacement = ',' + s
        protein = protein.replace(s, replacement)

    for s in sites:
        for n in no:
            anti_cleavage = s + ',' + n
            if pos == 'n':
                anti_cleavage = ',' + s + n
            protein = protein.replace(anti_cleavage, s + n)

    peptides = list(filter(lambda x: len(x) >= min_len, protein.split(',')))
    return [peptide for peptide in peptides if peptide]


def digest(protein, sites='KR', pos='c', no='P', min_len=0):
    if _FAST_DIGEST_ENABLED and _FAST_DIGEST_AVAILABLE:
        return _FAST_DIGEST_MODULE.digest(protein, sites, pos, no, min_len)
    return _digest_python(protein, sites, pos, no, min_len)


def _trypsin_python(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40):
    """
    Digest a protein sequence and return peptides within the requested length range.
    """
    peptides_cut_all = _digest_python(protein=protein, sites=sites, pos=pos, no=no, min_len=0)
    peptides = []
    for i in range(miss_cleavage + 1):
        for j in range(len(peptides_cut_all) - i):
            peptide = ''.join(peptides_cut_all[j:j + i + 1])
            peplen = len(peptide)
            if peplen >= peplen_min and peplen <= peplen_max:
                peptides.append(peptide)
    return peptides


def TRYPSIN(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40):
    if _FAST_DIGEST_ENABLED and _FAST_DIGEST_AVAILABLE:
        return _FAST_DIGEST_MODULE.trypsin(protein, sites, pos, no, miss_cleavage, peplen_min, peplen_max)
    return _trypsin_python(protein, sites, pos, no, miss_cleavage, peplen_min, peplen_max)


@lru_cache(maxsize=200000)
def get_new_pep_after_checkSimilar(seq, checkSimilar='GG=N,N=D,Q=E'):
    """
    Replace amino acids as described in checkSimilar from left to right.
    """
    seqnew = seq
    for replacement in checkSimilar.split(','):
        old, new = replacement.split('=')
        seqnew = seqnew.replace(old, new)
    return seqnew


def splitStringWithPeptide(proseq, peptide, anti_cleavage_sites='P', cleavage_sites='KR', pos='c'):
    """Split a protein sequence into segments using a peptide boundary."""
    if not peptide or peptide not in proseq:
        return [proseq]

    pep_from_pep = digest(protein=peptide, sites=cleavage_sites, pos=pos, no=anti_cleavage_sites, min_len=0)
    pep_from_pep = [item for item in pep_from_pep if item]
    n_miss_cleavage_site = len(pep_from_pep) - 1
    if n_miss_cleavage_site > 10:
        print(proseq, peptide, 'double check cleavage')

    peptides = digest(protein=proseq, sites=cleavage_sites, pos=pos, no=anti_cleavage_sites, min_len=0)
    peptides = [item for item in peptides if item]
    positions = [[0, len(peptides[0])]]
    for peptide_part in peptides[1:]:
        positions.append([positions[-1][1], positions[-1][1] + len(peptide_part)])

    positions_peptide = []
    for i in range(len(peptides) - n_miss_cleavage_site):
        pep_miss = ''.join(peptides[i:i + n_miss_cleavage_site + 1])
        pep_miss_pos = [positions[i][0], positions[i + n_miss_cleavage_site][1]]
        if pep_miss == peptide:
            positions_peptide.append(pep_miss_pos)

    while True:
        for i in range(len(positions_peptide) - 1):
            if positions_peptide[i][1] > positions_peptide[i + 1][0]:
                positions_peptide.pop(i + 1)
                break
        else:
            break

    break_points = [0] + [point for pair in positions_peptide for point in pair] + [len(proseq)]
    positions_peptide = [[break_points[i], break_points[i + 1]] for i in range(len(break_points) - 1)]
    peptides_seg = [proseq[i:j] for i, j in positions_peptide]
    return [item for item in peptides_seg if item]


def get_target_peptides(args):
    """Digest proteins from input files. Used when checkSimilar is enabled."""
    start_time = time.time()
    target_peptide_extra = set()
    for file_fasta in args.fasta:
        for _header, seq in read_fasta_file(file_path=file_fasta):
            target_protein_seq = seq
            if not args.iso:
                target_protein_seq = target_protein_seq.replace('I', 'L')
            target_peptide_extra.update(
                TRYPSIN(
                    target_protein_seq,
                    miss_cleavage=args.miss_cleavage,
                    peplen_min=args.minlen,
                    peplen_max=args.maxlen,
                    sites=args.csites,
                    no=args.noc,
                    pos=args.cpos,
                )
            )

    upeps_extra = target_peptide_extra
    print('target peptides after allowing missed clevage:', len(upeps_extra))
    print('checkSimilar is enabled!', args.checkSimilar)
    upeps_extra2 = set([get_new_pep_after_checkSimilar(i, args.checkSimilar) for i in upeps_extra])
    print('target peptides after', args.checkSimilar, len(upeps_extra2))
    print("--- {:.2f} seconds digest input proteins allowing miss cleavage ---".format(time.time() - start_time))
    args.upeps_extra2 = upeps_extra2


def get_decoy_peptides(args):
    """Digest proteins from the final decoy file. Used when checkSimilar is enabled."""
    start_time = time.time()
    decoy_peptide_extra = set()
    for _header, seq in read_fasta_file(file_path=args.dout):
        decoy_protein_seq = seq
        if not args.iso:
            decoy_protein_seq = decoy_protein_seq.replace('I', 'L')
        decoy_peptide_extra.update(
            TRYPSIN(
                decoy_protein_seq,
                miss_cleavage=args.miss_cleavage,
                peplen_min=args.minlen,
                peplen_max=args.maxlen,
                sites=args.csites,
                no=args.noc,
                pos=args.cpos,
            )
        )

    dpeps_extra = decoy_peptide_extra
    print('decoy peptides after allowing missed clevage:', len(dpeps_extra))
    print('checkSimilar is enabled!', args.checkSimilar)
    dpeps_extra2 = set([get_new_pep_after_checkSimilar(i, args.checkSimilar) for i in dpeps_extra])
    print('decoy peptides after ', args.checkSimilar, len(dpeps_extra2))
    print("--- {:.2f} seconds digest decoy proteins allowing miss cleavage ---".format(time.time() - start_time))
    args.dpeps_extra2 = dpeps_extra2


def build_parser():
    parser = argparse.ArgumentParser(
        description='Digest protein sequences from FASTA input or a direct sequence string.'
    )
    parser.add_argument(
        'fasta',
        metavar='*.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz',
        nargs='*',
        help='FASTA file(s) containing protein sequences to digest',
    )
    parser.add_argument(
        '--sequence',
        default='',
        help='Digest a single protein sequence directly instead of reading FASTA input',
    )
    parser.add_argument(
        '--cleavage_sites',
        '-c',
        dest='csites',
        default='KR',
        help='A list of amino acids at which to cleave during digestion. Default = KR',
    )
    parser.add_argument(
        '--anti_cleavage_sites',
        '-a',
        dest='noc',
        default='P',
        help='A list of amino acids at which not to cleave if following cleavage site. Default = P',
    )
    parser.add_argument(
        '--cleavage_position',
        '-p',
        dest='cpos',
        default='c',
        choices=['c', 'n'],
        help='Set cleavage to be c or n terminal of specified cleavage sites. Default = c',
    )
    parser.add_argument(
        '--min_peptide_length',
        '-l',
        dest='minlen',
        default=6,
        type=int,
        help='Set minimum peptide length. Default = 6',
    )
    parser.add_argument(
        '--max_peptide_length',
        '-M',
        dest='maxlen',
        default=40,
        type=int,
        help='Set maximum peptide length. Only used by trypsin mode. Default = 40',
    )
    parser.add_argument(
        '--miss_cleavage',
        '-L',
        dest='miss_cleavage',
        default=2,
        type=int,
        help='Missed cleavage count for trypsin mode. Default = 2',
    )
    parser.add_argument(
        '--no_isobaric',
        '-i',
        dest='iso',
        default=False,
        action='store_true',
        help='Do not change I to L before digestion. Default=False',
    )
    parser.add_argument(
        '--fast_digest',
        dest='fast_digest',
        default=False,
        action='store_true',
        help='Use Cython-accelerated digestion if available. Default=False',
    )
    parser.add_argument(
        '--method',
        choices=['digest', 'trypsin'],
        default='digest',
        help='Digestion method to run. Default=digest',
    )
    parser.add_argument(
        '--output-format',
        choices=['tsv', 'peptide', 'fasta'],
        default='tsv',
        help='Output format. Default=tsv',
    )
    parser.add_argument(
        '--header-template',
        default='{protein_id}_{index}',
        help='Header template for FASTA output. Default={protein_id}_{index}',
    )
    return parser


def _normalize_sequence(sequence, isobaric):
    sequence = sequence.upper().strip('*')
    if not isobaric:
        sequence = sequence.replace('I', 'L')
    return sequence


def _digest_sequence(sequence, args):
    if args.method == 'trypsin':
        return TRYPSIN(
            sequence,
            sites=args.csites,
            pos=args.cpos,
            no=args.noc,
            miss_cleavage=args.miss_cleavage,
            peplen_min=args.minlen,
            peplen_max=args.maxlen,
        )
    return digest(sequence, sites=args.csites, pos=args.cpos, no=args.noc, min_len=args.minlen)


def iter_digestion_results(args):
    if args.sequence:
        sequence = _normalize_sequence(args.sequence, args.iso)
        yield 'sequence', _digest_sequence(sequence, args)
        return

    for file_fasta in args.fasta:
        for header, sequence in read_fasta_file(file_path=file_fasta):
            normalized_sequence = _normalize_sequence(sequence, args.iso)
            yield header, _digest_sequence(normalized_sequence, args)


def _format_fasta_header(header, index, template):
    header_text = header[1:] if header.startswith('>') else header
    try:
        return template.format(protein_id=header_text, index=index, header=header_text)
    except KeyError as exc:
        raise ValueError('unsupported header-template field: {}'.format(exc.args[0]))


def write_output(results, args, handle):
    for header, peptides in results:
        if args.output_format == 'tsv':
            for peptide in peptides:
                handle.write('{}\t{}\n'.format(header, peptide))
            continue
        if args.output_format == 'peptide':
            for peptide in peptides:
                handle.write(peptide + '\n')
            continue
        for index, peptide in enumerate(peptides, start=1):
            handle.write('>{}\n{}\n'.format(_format_fasta_header(header, index, args.header_template), peptide))


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.sequence) == bool(args.fasta):
        parser.error('provide either FASTA input or --sequence')
    set_fast_digest(args.fast_digest)
    try:
        write_output(iter_digestion_results(args), args, sys.stdout)
    except ValueError as exc:
        parser.error(str(exc))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
