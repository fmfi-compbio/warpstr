import os
import re
from dataclasses import dataclass
from multiprocessing import Pool
from typing import List, Tuple

import numpy as np

import src.dtw_automata.StateAutomata as sta
import src.templates as tmpl
from src.dtw_automata.overview import load_overview, store_collapsed, store_results
from src.dtw_automata.plotter import plot_collapsed, plot_summaries
from src.extractor.tr_extractor import Flanks, load_flanks
from src.input_handler.fast5 import Fast5
from src.input_handler.input import caller_config, main_config
from src.input_handler.locus import Locus
from src.squiggler.pore_model import pore_model

from .caller import CallerResult, WarpSTR


def main_wrapper(locus: Locus):
    """
    Wrapper function for processing workload of all reads
    """
    overview_path, df_overview = load_overview(locus.path)
    workload = get_workload(df_overview, locus.path)

    dtw_automat = CallerWrapper(locus)
    results = dtw_automat.run_automata(workload)

    print(results)
    seq_results = [(i.seq, i.resc_seq) for i in results]
    cost_results = [(i.cost, i.resc_cost) for i in results]

    df_overview = store_results(overview_path, df_overview, seq_results, cost_results, locus.path)
    plot_summaries(locus.path, df_overview)

    if len(dtw_automat.units) > 0:
        collapsed_results = [dtw_automat.collapse_repeats(i[1]) for i in seq_results]
        reverse_lst = [i.reverse for i in workload]
        df_collapsed = store_collapsed(
            collapsed_results, dtw_automat.units, dtw_automat.repeat_units, reverse_lst, locus)
        plot_collapsed(df_collapsed, locus.path)

    return df_overview, df_collapsed


@dataclass
class ReadSignal:
    name: str
    reverse: bool
    signal: np.ndarray


def get_workload(df_overview, path: str) -> List[ReadSignal]:
    """
    Prepares all the data for processing
    """
    tr_regions: List[ReadSignal] = []
    for row in df_overview.itertuples():
        if row.saved:
            fast5path = os.path.join(path, tmpl.FAST5_SUBDIR, str(
                row.run_id), tmpl.ANNOT_SUBDIR, row.Index+'.fast5')
            fast5 = Fast5(fast5path, get_data_only=True)
            norm_signal = fast5.get_data_processed((row.l_start_raw, row.r_end_raw))
            tr_regions.append(ReadSignal(row.Index, row.reverse, norm_signal))
    return tr_regions


class CallerWrapper:
    temp_sta: sta.StateAutomata
    rev_sta: sta.StateAutomata

    def __init__(self, locus: Locus):
        self.locus = locus
        self.template_seq, self.reverse_seq = get_seqs(locus.sequence, load_flanks(locus.path))
        self.problems = self.check_high_similarity(locus.sequence)
        self.units, self.repeat_units, self.offsets = self.break_into_units(locus.sequence)
        print('units:', self.units)
        print('repeat_units:', self.repeat_units)
        print('offsets:', self.offsets)

        self.temp_sta = sta.StateAutomata(self.template_seq)
        self.rev_sta = sta.StateAutomata(self.reverse_seq)

    def init_pool(self):
        global out_warp_path
        global rev_sta
        global temp_sta
        global flank_length

        if caller_config.visualize_alignment:
            out_warp_path = os.path.join(self.locus.path, tmpl.PREDICTIONS_SUBDIR, tmpl.WARPS)
        else:
            out_warp_path = None
        rev_sta = self.rev_sta
        temp_sta = self.temp_sta
        flank_length = self.locus.flank_length

    def run_automata(self, workload: List[ReadSignal]):
        """
        Run WarpSTR automata for each piece of signal in workload
        """
        results = []
        print('running automata')
        if main_config.threads > 1:
            with Pool(main_config.threads, initializer=self.init_pool) as pool:
                results = pool.map(warpstr_call_parallel, workload)
        else:
            results = [warpstr_call_sequential(i, self.locus.flank_length) for i in workload]
        return results

    def check_high_similarity(self, sequence: str):
        """
        Checks high similarity of expected signal values
        """
        template = sequence
        reverse = reverse_uniq_sequence(sequence)

        diffs_t = self.get_diffs_for_all(template)
        diffs_r = self.get_diffs_for_all(reverse)
        print('diffs', diffs_t)
        out_path = os.path.join(self.locus.path, tmpl.SUMMARY_SUBDIR, 'state_similarity.csv')
        with open(out_path, 'w') as file:
            file.write('pattern,strand,mean_diff,median_diff\n')
            for i in diffs_t:
                file.write('{pattern},template,{m_diff:.3f},{med_diff:.3f}\n'.format(pattern=i,
                                                                                     m_diff=diffs_t[i][0], med_diff=diffs_t[i][1]))
            for i in diffs_r:
                file.write('{pattern},reverse,{m_diff:.3f},{med_diff:.3f}\n'.format(pattern=i,
                                                                                    m_diff=diffs_r[i][0], med_diff=diffs_r[i][1]))

        template_problems = []
        for i in diffs_t:
            if caller_config.min_state_similarity > diffs_t[i][0] or caller_config.min_state_similarity > diffs_t[i][1]:
                problem = {}
                problem['pattern'] = i
                problem['mean_diff'] = diffs_t[i][0]
                problem['median_diff'] = diffs_t[i][1]
                template_problems.append(problem)

        reverse_problems = []
        for i in diffs_r:
            if caller_config.min_state_similarity > diffs_r[i][0] or caller_config.min_state_similarity > diffs_r[i][1]:
                problem = {}
                problem['pattern'] = i
                problem['mean_diff'] = diffs_r[i][0]
                problem['median_diff'] = diffs_r[i][1]
                reverse_problems.append(problem)

        for i in template_problems:
            print('Warning: Template has repeat unit {} with high state similarity'.format(i['pattern']))
        for i in reverse_problems:
            print('Warning: high similarity of state values in reverse pattern {}'.format(i['pattern']))

        return (template_problems, reverse_problems)

    def break_into_units(self, template: str):
        """
        Breaks input config sequence into repeat units
        """
        que: List[int] = []
        units: List[str] = []
        offsets: List[int] = []
        offset = 0
        for idx, i in enumerate(template):
            if i in ('(', '{'):
                que.append(idx)
            elif i in (')', '}'):
                val = que.pop()
                if len(que) == 0:
                    units.append(template[val:idx+1])  # .strip("()\{\}"))
                    offsets.append(offset)
                    offset = 0
            elif len(que) == 0:
                offset += 1

        repeat_units: List[List[str]] = []
        for unit in units:
            repeats: List[str] = []
            seq = ''
            base_pattern = ''.join([c for c in unit if c not in ['(', ')']])
            for idx, i in enumerate(base_pattern):
                if i == '{':
                    repeats.append(seq)
                    seq = ''
                elif i == '}':
                    nl: List[str] = []
                    for p in repeats:
                        nl.append(p+seq)
                    for p in nl:
                        repeats.append(p)
                    seq = ''
                else:
                    seq += i
            if len(seq) > 0:
                repeats.append(seq)

            unwrapped: List[str] = []
            for rep in repeats:
                pattern = ['']
                for char in rep:
                    if char not in tmpl.DNA_DICT:
                        for idx, p in enumerate(pattern):
                            pattern[idx] = p + char
                    else:
                        new: List[str] = []
                        for iupac in tmpl.DNA_DICT[char]:
                            for p in pattern:
                                new.append(p+iupac)
                        pattern = new
                for p in pattern:
                    unwrapped.append(p)
            repeat_units.append(unwrapped)

        return units, repeat_units, offsets

    def get_diffs_for_all(self, template: str):
        """
        Gets differences between expected signals for each unit
        """
        diffs = {}
        brackets = ['(', ')', '{', '}']
        res = re.findall(r'[\(\{].*?[\)\}]', template)

        for r in res:
            base_pattern = ''.join([c for c in r if c not in brackets])
            pattern = ['']
            for char in base_pattern:
                if char not in tmpl.DNA_DICT:
                    for idx, p in enumerate(pattern):
                        pattern[idx] = p + char
                else:
                    new: List[str] = []
                    for iupac in tmpl.DNA_DICT[char]:
                        for p in pattern:
                            new.append(p+iupac)
                    pattern = new
            for p in pattern:
                diffs[p] = self.get_consecutive_diff(p)

        return diffs

    def get_consecutive_diff(self, pattern: str) -> Tuple[float, float]:
        """
        Get mean and median differences between expected signals
        """
        rep = pattern*pore_model.kmersize
        kmers = [rep[i:pore_model.kmersize+i] for i in range(len(pattern)+1)]
        pore_values = np.abs(np.diff(
            [pore_model.get_value(kmer) for kmer in kmers]))
        return np.mean(pore_values), np.median(pore_values)

    def collapse_repeats(self, seq: str):
        """
        Collapses repeat units
        """
        results = []
        slide = seq
        for i in self.repeat_units:
            if isinstance(i, list):
                results.append([0]*len(i))
            else:
                results.append(0)
        for idx, (p, off) in enumerate(zip(self.repeat_units, self.offsets)):
            slide = slide[off:]
            while len(slide) > 0:
                if (isinstance(p, str) and slide[:len(p)] == p):
                    results[idx] += 1
                    slide = slide[len(p):]
                elif isinstance(p, list):
                    notfound = True
                    for idx2, possible in enumerate(p):
                        if possible == slide[:len(possible)]:
                            results[idx][idx2] += 1
                            rep = slide[len(possible):]
                            notfound = False
                else:
                    break

                if notfound:
                    break
                slide = rep
        return results


# states, endstate, mask
def warpstr_call_sequential(input_data: ReadSignal, flank_length: int) -> CallerResult:
    # warper = Warper(flank_length, input_data.signal, states, endstate)
    return CallerResult(
        seq='TEST',
        resc_seq='TEST',
        cost=1.1,
        resc_cost=1.1
    )


def warpstr_call_parallel(input_data: ReadSignal):
    """
    :param input_data: list of lists with signal and reverse info
    :param rev_sta: global StateAutomata
    :param temp_sta: global StateAutomata
    :param flank_length: global int
    :param out_warp_path: global str
    :returns list of sequence lengths
    """

    if input_data.reverse:
        states = rev_sta.states
        endstate = rev_sta.endstate
        repeat_mask = rev_sta.mask
    else:
        states = temp_sta.states
        endstate = temp_sta.endstate
        repeat_mask = temp_sta.mask

    warpstr = WarpSTR(flank_length, states, endstate, repeat_mask, out_warp_path, input_data.reverse, input_data.name)
    return warpstr.run(input_data.signal)


def get_seqs(sequence: str, flanks: Flanks):
    """ Gets whole input sequence for state automaton using repeating part and flanks
    """
    tmp = flanks.template.left + sequence + flanks.template.right
    rev = flanks.reverse.left + reverse_uniq_sequence(sequence) + flanks.reverse.right
    return tmp, rev


def reverse_uniq_sequence(sequence: str):
    """ Transforms the regular expression to the reverse strand """
    rev = ''
    for i in sequence[::-1]:
        rev += tmpl.ENCODING_DICT[i]
    return rev
