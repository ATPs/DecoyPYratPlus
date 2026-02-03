# cython: language_level=3
# cython: boundscheck=False, wraparound=False, cdivision=True

from cpython.bytes cimport PyBytes_AS_STRING
from cpython.unicode cimport PyUnicode_AsASCIIString

cdef inline bytes _to_bytes(object s):
    if isinstance(s, bytes):
        return <bytes>s
    return (<bytes>PyUnicode_AsASCIIString(s))

cdef list _digest_spans(bytes seq, bytes sites, bytes no, char pos):
    cdef unsigned char cleave[256]
    cdef unsigned char nocut[256]
    cdef int i
    for i in range(256):
        cleave[i] = 0
        nocut[i] = 0
    cdef Py_ssize_t n = len(sites)
    for i in range(n):
        cleave[sites[i]] = 1
    n = len(no)
    for i in range(n):
        nocut[no[i]] = 1

    cdef Py_ssize_t seqlen = len(seq)
    cdef list spans = []
    cdef Py_ssize_t start = 0
    cdef Py_ssize_t cut
    cdef unsigned char c
    for i in range(seqlen):
        c = seq[i]
        if cleave[c] and (i + 1 == seqlen or not nocut[seq[i + 1]]):
            if pos == <char>99:  # ord('c')
                cut = i + 1
            else:
                cut = i
            if cut > start:
                spans.append((start, cut))
            start = cut
    if start < seqlen:
        spans.append((start, seqlen))
    return spans


def digest(protein, sites='KR', pos='c', no='P', min_len=0):
    cdef bytes seq = _to_bytes(protein)
    cdef bytes sites_b = _to_bytes(sites)
    cdef bytes no_b = _to_bytes(no)
    cdef bytes pos_b = _to_bytes(pos)
    cdef char pos_c = pos_b[0]
    cdef list spans = _digest_spans(seq, sites_b, no_b, pos_c)
    cdef list peptides = []
    cdef Py_ssize_t i
    cdef Py_ssize_t start
    cdef Py_ssize_t end
    for i in range(len(spans)):
        start = spans[i][0]
        end = spans[i][1]
        if end - start >= min_len:
            peptides.append(seq[start:end].decode('ascii'))
    return peptides


def trypsin(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40):
    cdef bytes seq = _to_bytes(protein)
    cdef bytes sites_b = _to_bytes(sites)
    cdef bytes no_b = _to_bytes(no)
    cdef bytes pos_b = _to_bytes(pos)
    cdef char pos_c = pos_b[0]
    cdef list spans = _digest_spans(seq, sites_b, no_b, pos_c)
    cdef list peptides = []
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef Py_ssize_t start
    cdef Py_ssize_t end
    cdef Py_ssize_t length
    cdef Py_ssize_t nspans = len(spans)
    for i in range(miss_cleavage + 1):
        for j in range(nspans - i):
            start = spans[j][0]
            end = spans[j + i][1]
            length = end - start
            if length >= peplen_min and length <= peplen_max:
                peptides.append(seq[start:end].decode('ascii'))
    return peptides
