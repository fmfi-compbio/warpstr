import os
from typing import List

import numpy as np
import pandas as pd

import src.templates as tmpl
from src.input_handler.locus import Locus


def store_collapsed(results, units: List[str], rep_units: List[List[str]], reverse_lst: List[bool], locus: Locus):
    """
    Stores results as given by repeat units
    """
    preds = {}
    for idx, i in enumerate(units):
        if len(results[0][idx]) > 1:
            main = 'main_'+rep_units[idx][0]
            preds[main] = np.array([np.sum(j[idx]) for j in results])
            for idx2, k in enumerate(rep_units[idx][1:]):
                inter = 'inter_'+k[len(rep_units[idx][0]):]
                preds[inter] = np.array([j[idx][idx2+1] for j in results])
        else:
            name = i.strip('(').strip(')')
            preds[name] = np.array([j[idx][0] for j in results])
    preds['reverse'] = reverse_lst
    df_preds = pd.DataFrame.from_dict(preds)

    out_path = os.path.join(
        locus.path, tmpl.PREDICTIONS_SUBDIR, tmpl.COMPLEX_SUBDIR)
    if os.path.isdir(out_path) is False:
        os.mkdir(out_path)
    df_preds.to_csv(os.path.join(out_path, 'complex_repeat_units.csv'))
    return df_preds


def load_overview(locus_path: str):
    """
    Loads overview file as dataframe
    """
    overview_path = os.path.join(locus_path, tmpl.OVERVIEW_NAME)
    try:
        df_overview = pd.read_csv(overview_path)
        df_overview.set_index('read_name', inplace=True)
        df_overview.columns = df_overview.columns.map(str)
    except FileNotFoundError:
        raise FileNotFoundError(f'Not found the overview file {overview_path} - Please check the "output" in config')
    return overview_path, df_overview


def store_results(overview_path, df_overview, seq_results, cost_results, locus_path: str):
    """
    Stores results in overview file
    """
    fasta_lst, newcol, dbg1, dbg2, dbg3 = append_results(
        seq_results, cost_results, df_overview)
    write_results_to_fasta(fasta_lst, locus_path)
    df_overview = save_overview(
        overview_path, df_overview, newcol, dbg1, dbg2, dbg3)
    return df_overview


def append_results(seq_results, cost_results, df_overview):
    """
    Appends results
    """
    fasta_lst, newcol, dbg1, dbg2, dbg3 = [], [], [], [], []
    idx = 0
    for row in df_overview.itertuples():
        if row.saved:
            newcol.append(len(seq_results[idx][1]))
            dbg1.append(len(seq_results[idx][0]))
            dbg2.append(cost_results[idx][0])
            dbg3.append(cost_results[idx][1])
            fasta_lst.append((row.Index, seq_results[idx][1], row.reverse))
            idx = idx + 1
        else:
            newcol.append(-1)
            dbg1.append(-1)
            dbg2.append(-1)
            dbg3.append(-1)
    return fasta_lst, newcol, dbg1, dbg2, dbg3


def write_results_to_fasta(fasta_lst, locus_path: str):
    """
    Stores results as FASTA sequences
    """
    fasta_out = os.path.join(
        locus_path, tmpl.PREDICTIONS_SUBDIR, 'sequences', 'all.fasta')
    fasta_out_template = os.path.join(
        locus_path, tmpl.PREDICTIONS_SUBDIR, 'sequences', 'sequences_template.fasta')
    fasta_out_reverse = os.path.join(
        locus_path, tmpl.PREDICTIONS_SUBDIR, 'sequences', 'sequences_reverse.fasta')
    with open(fasta_out, 'w') as file:
        for fastaid, fastaseq, rev in fasta_lst:
            file.write('>'+fastaid+'\n')
            file.write(fastaseq+'\n')
            file.write('\n')
    with open(fasta_out_template, 'w') as file:
        for fastaid, fastaseq, rev in fasta_lst:
            if rev is False:
                file.write('>'+fastaid+'\n')
                file.write(fastaseq+'\n')
                file.write('\n')
    with open(fasta_out_reverse, 'w') as file:
        for fastaid, fastaseq, rev in fasta_lst:
            if rev:
                file.write('>'+fastaid+'\n')
                file.write(fastaseq+'\n')
                file.write('\n')


def save_overview(overview_path, df_overview, newcol, dbg1, dbg2, dbg3):
    """
    Save Warping results in overview
    """
    prev_res = [c for c in df_overview.columns if c.startswith('result')]
    df_overview.drop(columns=prev_res, inplace=True)

    df_overview['results'] = newcol
    df_overview['orig'] = dbg1
    df_overview['dtw_cost1'] = dbg2
    df_overview['dtw_cost2'] = dbg3

    df_overview.to_csv(overview_path)
    print(f' Results stored in overview file {overview_path}')
    return df_overview
