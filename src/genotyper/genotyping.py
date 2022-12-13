import os
from typing import List, Optional

import numpy as np
import pandas as pd
from sklearn.mixture import BayesianGaussianMixture

import src.templates as tmpl
from src.caller.overview import load_overview
from src.config import genotyping_config
from src.schemas import Genotype, find_nearest

from .muscle import run_muscle
from .plotter import plot_clustering_preds, plot_complex_repeats


def decode_alleles_complex(gmm_out_dict, df):
    if gmm_out_dict['is_hetero']:
        df1 = df.loc[gmm_out_dict['group1']]
        df2 = df.loc[gmm_out_dict['group2']]
        mediangroup1 = []
        for col in df1.columns:
            if col != 'reverse':
                mediangroup1.append(find_nearest(df1[col], np.median(df1[col])))
        mediangroup1_cnt = len(gmm_out_dict['group1'])
        mediangroup2 = []
        for col in df2.columns:
            if col != 'reverse':
                mediangroup2.append(find_nearest(df2[col], np.median(df2[col])))
        mediangroup2_cnt = len(gmm_out_dict['group2'])
    else:
        df1 = df.loc[gmm_out_dict['group1']]
        df2 = None
        mediangroup1 = []
        for col in df1.columns:
            if col != 'reverse':
                mediangroup1.append(find_nearest(df1[col], np.median(df1[col])))
        mediangroup1_cnt = len(gmm_out_dict['group1'])
        mediangroup2 = '-'
        mediangroup2_cnt = '-'
    return (mediangroup1, mediangroup1_cnt, mediangroup2, mediangroup2_cnt)


def filter_out_complex(df, cols, STD_COEFF: float):
    if len(df) <= 5:
        return df

    idxes = []
    for col in cols:
        mean = np.mean(df[col])
        std = np.std(df[col])
        idxes.append([idx for idx, i in enumerate(df[col]) if i < (mean-STD_COEFF*std) or i > (mean+STD_COEFF*std)])

    idxes = np.unique([item for sublist in idxes for item in sublist])

    df.drop(idxes, inplace=True)
    return df


def filter_out(values: List[int], STD_COEFF: float):
    if len(values) <= 5:
        return values

    mean = np.mean(values)
    std = np.std(values)
    return [i for i in values if i >= (mean-STD_COEFF*std) and i <= (mean+STD_COEFF*std)]


def run_genotyping_overview(overview, locus_path: str, muscle_path: Optional[str]):
    if overview is None:
        overview_path, overview = load_overview(locus_path)

    unfilt_values, unfilt_basecalls, use_basecalls = load_predictions(overview)

    gt = run_genotyping(unfilt_values)
    gt_bc = run_genotyping(unfilt_basecalls) if use_basecalls else None
    store_predictions(gt, gt_bc, locus_path)

    if genotyping_config.visualize:
        visualize_predictions(locus_path, gt, gt_bc)

    if muscle_path and genotyping_config.msa:
        run_muscle(locus_path, muscle_path, gt, gt_bc)


def visualize_predictions(locus_path: str, gt: Genotype, gt_bc: Optional[Genotype]):
    img_path = os.path.join(locus_path, tmpl.SUMMARY_SUBDIR, 'alleles.svg')
    vals = (gt.group1, gt.group2)
    if gt_bc:
        bvals = (gt_bc.group1, gt_bc.group2)
        plot_clustering_preds(img_path, vals, bvals, gt.alleles, gt_bc.alleles)
    else:
        plot_clustering_preds(img_path, vals, None, gt.alleles, None)


def load_predictions(overview):
    unfilt_values, unfilt_basecalls = [], []
    use_basecalls = ('r_seq_start' in overview.columns and 'l_seq_end' in overview.columns)
    for row in overview.itertuples():
        if row.saved:
            unfilt_values.append(row.results)
            if use_basecalls:
                unfilt_basecalls.append(row.r_seq_start-row.l_seq_end)
    return unfilt_values, unfilt_basecalls, use_basecalls


def store_predictions(gt: Genotype, gt_bc: Optional[Genotype], locus_path: str):
    final_preds_file = os.path.join(locus_path, 'predictions', 'alleles.csv')
    with open(final_preds_file, 'w') as f:
        f.write('WarpSTR_allele1,WarpSTR_allele1_freq,WarpSTR_allele2,WarpSTR_allele2_freq,'
                'basecall_allele1,basecall_allele1_freq,'
                'basecall_allele2,basecall_allele2_freq\n')
        f.write(f'{gt.first_allele},{gt.first_allele_sz},{gt.second_allele},{gt.second_allele_sz},')
        if gt_bc:
            f.write(f'{gt_bc.first_allele},{gt_bc.first_allele_sz},{gt_bc.second_allele},{gt_bc.second_allele_sz}')

    print(f'Allele lengths as given by WarpSTR: {gt.alleles}')
    if gt_bc:
        print(f'Allele lengths as given by basecall: {gt_bc.alleles}')


def run_genotyping_complex(locus_path: str, df):
    if df is None:
        inpath = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR, tmpl.COMPLEX_SUBDIR, 'complex_repeat_units.csv')
        if os.path.isfile(inpath):
            df = pd.read_csv(inpath, index_col=0)
        else:
            return

    cols = [col for col in df.columns if col != 'reverse']
    if len(cols) < 2:
        return

    df = filter_out_complex(df, cols, genotyping_config.std_filter)
    X = np.array(df[cols]).reshape(-1, len(cols))
    model_bayes = run_bayes(X)

    out = {}
    if one_large_group(model_bayes.weights_):
        out['is_hetero'] = False
        out['group1'] = df.index
        out['group2'] = []
        out['predictions'] = None
    else:
        preds = model_bayes.predict(X)
        out['is_hetero'] = True
        out['group1'] = [idx for idx, g in zip(df.index, preds) if g == 0]
        out['group2'] = [idx for idx, g in zip(df.index, preds) if g == 1]
        out['predictions'] = preds

    alleles = decode_alleles_complex(out, df)
    if out['predictions'] is not None:
        df['allele'] = preds

    img_path = os.path.join(locus_path, tmpl.SUMMARY_SUBDIR, 'complex_genotypes.svg')
    plot_complex_repeats(df, cols, alleles, img_path)


def run_genotyping(unfilt_vals: List[int]):
    vals = filter_out(unfilt_vals, genotyping_config.std_filter)

    if len(np.unique(vals)) == 1:
        return Genotype(group1=vals)

    X = np.array(vals).reshape(-1, 1)
    model_bayes = run_bayes(X)

    if one_large_group(model_bayes.weights_):
        return Genotype(group1=vals)

    preds: List[int] = model_bayes.predict(X)
    return Genotype(
        group1=[i for i, g in zip(vals, preds) if g == 0],
        group2=[i for i, g in zip(vals, preds) if g == 1],
        predictions=preds
    )


def one_large_group(weights: List[float]) -> bool:
    return any(x < genotyping_config.min_weight for x in weights)


def run_bayes(X: np.ndarray):
    return BayesianGaussianMixture(
        weight_concentration_prior=0.25,
        covariance_type='tied',
        n_components=2,
        n_init=5,
        max_iter=1000
    ).fit(X)
