
from re import A
import pysam
import argparse
import pandas as pd
import yaml
from collections import OrderedDict

def prepare_sequence(locus, reference_path):
    fa_seq = pysam.faidx(reference_path,locus['coord'])
    seq = ''.join(fa_seq.split("\n")[1:])
    noting_seq = ''.join(fa_seq.split("\n")[1:])
    if len(seq)<=1:
        raise ValueError(f'Reference repeat sequence is empty for {locus["name"]} - check coord {locus["coord"]}')
    if 'motif' not in locus:
        raise ValueError(f'Motif is not defined for {locus["name"]}. Define either motif or sequence')

    for motif in locus['motif'].split(','):
        a = '-'.join(seq.split(motif))
        counter = 0
        out_sequence = ''
        for ch in a:
            if ch=='-':
                counter += 1
            else:
                if counter>1:
                    out_sequence += f'({motif})'
                elif counter==1:
                    out_sequence += f'{motif}'
                out_sequence += ch
                counter = 0
        if counter>0:
            if counter==1:
                out_sequence += f'{motif}'
            else:
                out_sequence += f'({motif})'
        
        counter = 0
        a = '-'.join(noting_seq.split(motif))
        noting = ''
        for ch in a:
            if ch=='-':
                counter += 1
            else:
                if counter>1:
                    noting += f'({motif})[{counter}]'
                elif counter==1:
                    noting += f'{motif}'
                noting += ch
                counter = 0
        if counter>0:
            if counter==1:
                noting += f'{motif}'
            else:
                noting += f'({motif})[{counter}]'

        seq = out_sequence
        noting_seq = noting
        
    locus['noting'] = noting_seq
    locus['sequence'] = out_sequence
    return locus