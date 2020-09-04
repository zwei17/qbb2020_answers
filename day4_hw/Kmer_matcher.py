#!/usr/bin/env python3


import sys
from fasta_parser import FASTAReader


def kmer_matcher(target, query, k, head = False, Num = 1):
    
    target = open(target)
    target.seek(0,0)
    query = open(query)
    query.seek(0, 0)
    
    t_seqs = FASTAReader(target)
    q_seq = FASTAReader(query)[0][1]
    kmer_matches = []


    k = int(k)
    Num = int(Num)
    q = 0
    n = 0
    for q in range(len(q_seq) - k + 1):
        q_kmer = q_seq[q:(q + k)].upper()
        total = (len(q_seq) - k + 1)
        for seq_id, t_each in t_seqs:
            i = 0
            for i in range(len(t_each) - k + 1):
                if t_each[i].upper() != q_kmer[0]:
                    continue
                if q_kmer == t_each[i:(i + k)].upper():
                    kmer_matches.append((seq_id, i, q))
                    n += 1
                    print('{}: {}---{}'.format(seq_id, i, q))
                    if head and n == Num:
                        break
                if t_each[i+k-1] != q_kmer[0]:
                    i += k-1
            if head and n == Num:
                        break
        if head and n == Num:
            break


    target.close()
    query.close()
    return kmer_matches

if __name__ == '__main__':
    seqs = kmer_matcher(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    for seq_id, i, q in seqs:
        print('{}: {}---{}'.format(seq_id, i, q), file = sys.stdout)