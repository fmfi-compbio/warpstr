import os
from typing import List, Optional, Tuple

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

import src.templates as tmpl
from src.schemas import Genotype


def run_muscle(locus_path: str, muscle_path: str, gt: Genotype, gt_bc: Optional[Genotype]):
    path_our = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR, 'sequences')
    trcalls = call_muscle(muscle_path, path_our)

    if gt.is_hetero:
        groups_from_pred(muscle_path, path_our, gt, trcalls)

    if gt_bc:
        path_bc = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR, 'basecalls')
        basecalls = call_muscle(muscle_path, path_bc)

        if gt_bc.is_hetero:
            groups_from_pred(muscle_path, path_bc, gt_bc, basecalls)


def call_muscle(muscle_path: str, path: str):
    if os.path.exists(muscle_path) is False:
        raise ValueError(f'Could not load the MUSCLE tool from {muscle_path}')
    calls_file = os.path.join(path, 'all.fasta')
    msa_out = os.path.join(path, 'msa_all.fasta')
    calls = [(str(record.id), str(record.seq)) for record in SeqIO.parse(calls_file, 'fasta')]
    muscle_cline = MuscleCommandline(muscle_path, input=calls_file, out=msa_out, verbose=False, quiet=True)
    os.system(str(muscle_cline))
    return calls


def groups_from_pred(muscle_path: str, path_out: str, gt: Genotype, calls: List[Tuple[str, str]]):
    groups = gt.predictions
    group1_seqs = [i for i, g in zip(calls, groups) if g == 0]
    group2_seqs = [i for i, g in zip(calls, groups) if g == 1]

    trcalls_file_group1 = os.path.join(path_out, 'group1.fasta')
    store_fasta(trcalls_file_group1, group1_seqs)

    msa_out_group1 = os.path.join(path_out, 'msa_group1.fasta')
    muscle_cline = MuscleCommandline(muscle_path, input=trcalls_file_group1,
                                     out=msa_out_group1, verbose=False, quiet=True)
    os.system(str(muscle_cline))

    trcalls_file_group2 = os.path.join(path_out, 'group2.fasta')
    store_fasta(trcalls_file_group2, group2_seqs)

    msa_out_group2 = os.path.join(path_out, 'msa_group2.fasta')
    muscle_cline = MuscleCommandline(muscle_path, input=trcalls_file_group2,
                                     out=msa_out_group2, verbose=False, quiet=True)
    os.system(str(muscle_cline))


def store_fasta(filename: str, fasta: List[Tuple[str, str]]):
    with open(filename, 'w') as file:
        for seqid, seqrec in fasta:
            file.write('>'+seqid+'\n')
            file.write(str(seqrec)+'\n\n')
