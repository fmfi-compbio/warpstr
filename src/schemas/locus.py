import os
from typing import Optional

import pysam

from src.config import main_config


class Locus:
    name: str
    coord: str
    sequence: str
    path: str
    flank_length: int
    motif: Optional[str] = None
    noting: Optional[str] = None

    def __init__(
        self,
        name: str,
        coord: str,
        sequence: Optional[str] = None,
        motif: Optional[str] = None,
        noting: Optional[str] = None,
        flank_length: Optional[int] = None

    ):
        self.name = name
        self.coord = coord
        self.path = os.path.join(main_config.output, self.name)
        self.motif = motif

        if not flank_length:
            self.flank_length = main_config.flank_length
            print(f'Flank length not set for locus - using default value of {main_config.flank_length}')
        else:
            self.flank_length = flank_length

        if sequence:
            self.sequence = sequence.upper()
            self.noting = noting

        else:
            seq, noting = self.prepare_sequence(main_config.reference_path)
            self.noting = noting
            self.sequence = seq.upper()

    def prepare_sequence(self, reference_path: str):
        fa_seq: str = pysam.faidx(reference_path, self.coord)
        seq = ''.join(fa_seq.split('\n')[1:])
        seq = seq.upper()
        noting_seq = seq
        if len(seq) <= 1:
            raise ValueError(f'Reference repeat sequence is empty for {self.name} - check coord {self.coord}')
        if not self.motif:
            raise ValueError(f'Motif is not defined for {self.name}. Define either motif or sequence')

        for motif in self.motif.split(','):
            a = '-'.join(seq.split(motif))
            counter = 0
            out_sequence = ''
            for ch in a:
                if ch == '-':
                    counter += 1
                else:
                    if counter > 1:
                        out_sequence += f'({motif})'
                    elif counter == 1:
                        out_sequence += f'{motif}'
                    out_sequence += ch
                    counter = 0
            if counter > 0:
                if counter == 1:
                    out_sequence += f'{motif}'
                else:
                    out_sequence += f'({motif})'

            counter = 0
            a = '-'.join(noting_seq.split(motif))
            noting = ''
            for ch in a:
                if ch == '-':
                    counter += 1
                else:
                    if counter > 1:
                        noting += f'({motif})[{counter}]'
                    elif counter == 1:
                        noting += f'{motif}'
                    noting += ch
                    counter = 0
            if counter > 0:
                if counter == 1:
                    noting += f'{motif}'
                else:
                    noting += f'({motif})[{counter}]'

            seq = out_sequence
            noting_seq = noting

        return seq, noting_seq
