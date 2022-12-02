import os
import re
from math import sqrt
from datetime import datetime
from multiprocessing import Pool
import numpy as np
from scipy import interpolate
from mappy import revcomp

import src.dtw_automata.StateAutomata as sta
from src.dtw_automata.plotter import save_warp_img, plot_summaries, plot_collapsed
from src.dtw_automata.overview import store_results, load_overview, store_collapsed
from src.squiggler.pore_model import load_pore_model
from src.extractor.tr_extractor import load_flanks
from src.input_handler.fast5 import Fast5
import src.templates as tmpl


def main_wrapper(locus_path, pore_model_path, locus, config, threads):
    """
    Wrapper function for processing workload of all reads
    """
    overview_path, df_overview = load_overview(locus_path)
    workload = get_workload(df_overview, locus_path, config)

    reverse_lst = [i[1] for i in workload]
    dtw_automat = DtwAutomata(locus_path, pore_model_path, locus, config, threads)
    results = dtw_automat.run_automata(workload)

    seq_results = [i[0] for i in results]
    cost_results = [i[1] for i in results]

    df_overview = store_results(locus_path, overview_path, df_overview, seq_results, cost_results)
    plot_summaries(locus_path, df_overview, config)

    if len(dtw_automat.units)>0:
        collapsed_results = [dtw_automat.collapse_repeats(i[1]) for i in seq_results]
        df_collapsed = store_collapsed(
            collapsed_results, dtw_automat.units, dtw_automat.repeat_units, reverse_lst, locus_path)
        plot_collapsed(df_collapsed,locus_path)

    return seq_results, df_overview, df_collapsed

class Warper:
    """
    Class for aligning signal with state automaton
    """
    def __init__(self,template,states,min_per_state,endstate,reverse, mask=None):
        self.template = template
        self.states = states
        self.min_per_state = min_per_state
        self.endstate = endstate
        self.reverse = reverse
        self.tlen, self.qlen = len(template), len(states)
        self.dtw_matrix = self.__init_matrix()
        self.trace = []
        self.seq = ""
        if mask is None:
            self.mask = [False]*len(template)
        else:
            self.mask = mask

        self.__calc_dtw_astates()
        self.__backtracking()
        self.cost = self.__get_norm_cost()
        self.state_transitions = self.trace[np.insert(np.diff(self.trace).astype(np.bool), 0, True)]
        self.warped_signal = self.__get_warped_signal()
        self.__extract_sequence()

    def __calc_cost(self,val1,val2):
        return abs(val1-val2)

    def __init_matrix(self):
        return np.full((self.tlen, self.qlen), np.inf)

    def __visited(self,value):
        return False if value == np.inf else True

    def __calc_dtw_astates(self):

        #calc first state
        start_state = self.states[0]
        start_val = self.__calc_cost(self.template[0],start_state.value)
        self.dtw_matrix[0,0] = start_val

        for i in range(1,self.min_per_state+1):
            self.dtw_matrix[0,i] = start_val+self.__calc_cost(self.template[i],start_state.value)

        #run dtw matrix
        boundary = FLANK_LENGTH-10
        after_repeat = self.states[-1].simple_idx-boundary
        first_threshold = 6*boundary
        second_threshold = len(self.template)-(6*boundary)
        for i,val in enumerate(self.template[self.min_per_state:],start=self.min_per_state):
            back = self.min_per_state - 1 if self.mask[i] else self.min_per_state
            for j,curr_st in enumerate(self.states):
                if i<first_threshold: #at start bound the alignment
                    if i<curr_st.simple_idx*4 and i>curr_st.simple_idx*15:
                        continue
                elif (i>second_threshold and curr_st.simple_idx<after_repeat):
                    continue

                if self.__visited(self.dtw_matrix[i-1,j]):
                    cost = self.dtw_matrix[i-1,j] + self.__calc_cost(val,curr_st.value)
                    if cost<self.dtw_matrix[i,j]:
                        self.dtw_matrix[i][j] = cost

                for prev in curr_st.incoming:
                    prev_idx = prev.realidx

                    if self.__visited(self.dtw_matrix[i-back,prev_idx]) is False:
                        continue
                    prev_st_val = prev.value
                    skip_cost = self.dtw_matrix[i-back,prev_idx]

                    for skip in self.template[i-back+1:i]:
                        skip_cost += self.__calc_cost(skip,prev_st_val)

                    skip_cost += self.__calc_cost(val,curr_st.value)
                    if skip_cost<self.dtw_matrix[i,j]:
                        self.dtw_matrix[i,j] = skip_cost

    def __backtracking(self):
        """
        Gets the best path through the states
        """
        newidx = self.endstate
        last = len(self.dtw_matrix)-1

        while last != 0:
            curr_dist = self.dtw_matrix[last][newidx]
            curr_state = self.states[newidx]

            if self.dtw_matrix[last-1][newidx] == np.inf:
                shift_delta = np.inf
            else:
                shift = self.dtw_matrix[last-1][newidx] + \
                    self.__calc_cost(self.template[last], curr_state.value)
                shift_delta = abs(shift - curr_dist)

            skip_delta = np.inf
            if self.mask[last]:
                back = self.min_per_state - 1
            else:
                back = self.min_per_state

            for prev in curr_state.incoming:
                prev_idx = prev.realidx

                #HANDLE COMING FROM NORMAL STATE
                if self.dtw_matrix[last-back][prev_idx] == np.inf:
                    continue
                skip_cost = self.dtw_matrix[last-back][prev_idx]
                for skip in self.template[last-back+1:last]:
                    skip_cost += self.__calc_cost(skip,prev.value)

                skip_cost += self.__calc_cost(self.template[last],curr_state.value)
                delta = abs(skip_cost-curr_dist)

                if delta<skip_delta:
                    skip_delta = delta
                    skip_idx = prev_idx

            if skip_delta<shift_delta:
                last = last - back
                self.trace.append(newidx)
                for j in range(back-1):
                    self.trace.append(skip_idx)
                newidx = skip_idx
            else:
                last = last - 1
                self.trace.append(newidx)

        self.trace.append(newidx)
        self.trace.reverse()
        self.trace = np.asarray(self.trace,dtype=int)

    def __get_norm_cost(self):
        """
        Calculate normalized DTW distance
        :return: normalized DTW distance
        """
        return self.dtw_matrix[-1][self.endstate]/len(self.trace)

    def __get_warped_signal(self):
        """
        Calculate warped signal using list of traversed k-mer states
        :return: list of pore model values for traversed k-mer states
        """
        return [self.states[i].value for i in self.trace]

    def __extract_sequence(self):
        """
        Calculates sequence derived from traversal of states
        :return: sequence given by WarpSTR
        """
        self.seq = ""

        for i in self.state_transitions:
            self.seq += self.states[i].kmer[-1]

        self.offset = self.states[self.state_transitions[0]].simple_idx
        self.seq = self.seq[FLANK_LENGTH-self.offset:-FLANK_LENGTH]

        if self.reverse:
            self.seq = revcomp(self.seq)

    def create_alignment(self, method, reps_as_one):
        """
        Creates alignment between signal values and signal states
        """
        if method == "mean":
            method = np.mean
        elif method == "median":
            method = np.median

        state_raws,state_val, state_exp = [] ,[], []
        if reps_as_one:
            state_indices = np.unique(self.state_transitions)
            for state_index in state_indices:
                raws = np.take(self.template,np.where(self.trace==state_index)[0])
                state_raws.append(raws)
                state_val.append(method(raws))
                state_exp.append(self.states[state_index].value)
            return state_raws,state_val, state_exp

        else:
            idx = 0
            for state_transition in self.state_transitions:
                raws = []
                for trace_elem in self.trace[idx:]:
                    if trace_elem!=state_transition:
                        break
                    raws.append(self.template[idx])
                    idx += 1
                state_raws.append(raws)
                state_val.append(method(raws))
                state_exp.append(self.states[state_transition].value)

            return state_raws,state_val, state_exp

def get_workload(df_overview, locus_path, config):
    """
    Prepares all the data for processing
    """
    tr_regions = []
    for row in df_overview.itertuples():
        if row.saved:
            fast5path = os.path.join(locus_path, tmpl.FAST5_SUBDIR, str(
                row.run_id), tmpl.ANNOT_SUBDIR, row.Index+".fast5")
            fast5 = Fast5(fast5path,config,get_data_only=True)
            norm_signal = fast5.get_data_processed((row.l_start_raw,row.r_end_raw))
            tr_regions.append((norm_signal, row.reverse, row.Index))
    return tr_regions

class DtwAutomata:
    """
    Class Automata
    """
    def __init__(self, locus_path, pore_model_path, locus, config, threads):
        self.config = config
        self.locus_path = locus_path
        self.flank_length = locus['flank_length']
        self.template_seq, self.reverse_seq = get_seqs(locus['sequence'], load_flanks(locus_path))
        self.pore_model, self.kmer_size = load_pore_model(pore_model_path)
        self.problems = self.check_high_similarity(locus['sequence'])
        self.units,self.repeat_units, self.offsets = self.break_into_units(locus['sequence'])
        self.temp_sta = sta.StateAutomata(self.pore_model, self.template_seq, self.kmer_size)
        self.rev_sta = sta.StateAutomata(self.pore_model, self.reverse_seq, self.kmer_size)
        self.threads = threads

    def run_automata(self, workload):
        """
        Run WarpSTR automata for each piece of signal in workload
        """
        global out_warp_path
        global min_values_per_state
        global rev_states
        global temp_states
        global rescaling_config
        global FLANK_LENGTH
        global temp_end_state
        global rev_end_state
        global rev_mask
        global temp_mask
        global tr_unit
        global visualize

        tr_unit = self.config['states_in_segment']
        out_warp_path = os.path.join(self.locus_path,tmpl.PREDICTIONS_SUBDIR,tmpl.WARPS)
        min_values_per_state = self.config['min_values_per_state']
        rev_states = self.rev_sta.kmer_states
        temp_states = self.temp_sta.kmer_states
        rev_mask = self.rev_sta.mask
        temp_mask = self.temp_sta.mask
        rescaling_config = self.config['rescaling']
        FLANK_LENGTH = self.flank_length
        temp_end_state = self.temp_sta.endstate
        rev_end_state = self.rev_sta.endstate
        visualize = self.config['visualize_alignment']

        results = []
        if self.threads>1:
            with Pool(self.threads) as pool:
                results = pool.map(dtw_wrapper,workload)
        else:
            results = [dtw_wrapper(i) for i in workload]
        return results

    def check_high_similarity(self, sequence):
        """
        Checks high similarity of expected signal values
        """
        template = sequence
        reverse = reverse_uniq_sequence(sequence)

        diffs_t = self.get_diffs_for_all(template)
        diffs_r = self.get_diffs_for_all(reverse)

        out_path = os.path.join(self.locus_path,tmpl.SUMMARY_SUBDIR,"state_similarity.csv")
        with open(out_path,"w") as file:
            file.write('pattern,strand,mean_diff,median_diff\n')
            for i in diffs_t:
                file.write("{pattern},template,{m_diff:.3f},{med_diff:.3f}\n".format(pattern=i,\
                    m_diff=diffs_t[i][0],med_diff=diffs_t[i][1]))
            for i in diffs_r:
                file.write("{pattern},reverse,{m_diff:.3f},{med_diff:.3f}\n".format(pattern=i,\
                    m_diff=diffs_r[i][0],med_diff=diffs_r[i][1]))

        template_problems = []
        for i in diffs_t:
            if self.config["min_state_similarity"]>diffs_t[i][0] or self.config["min_state_similarity"]>diffs_t[i][1]:
                problem = {}
                problem["pattern"] = i
                problem["mean_diff"] = diffs_t[i][0]
                problem["median_diff"] = diffs_t[i][1]
                template_problems.append(problem)

        reverse_problems = []
        for i in diffs_r:
            if self.config["min_state_similarity"]>diffs_r[i][0] or self.config["min_state_similarity"]>diffs_r[i][1]:
                problem = {}
                problem["pattern"] = i
                problem["mean_diff"] = diffs_r[i][0]
                problem["median_diff"] = diffs_r[i][1]
                reverse_problems.append(problem)

        for i in template_problems:
            print("Warning: Template has repeat unit {} with high state similarity".format(i["pattern"]))
        for i in reverse_problems:
            print("Warning: high similarity of state values in reverse pattern {}".format(i["pattern"]))

        return (template_problems, reverse_problems)

    def break_into_units(self,template):
        """
        Breaks input config sequence into repeat units
        """
        que = []
        units = []
        offsets = []
        offset = 0
        for idx,i in enumerate(template):
            if i in ('(', '{'):
                que.append(idx)
            elif i in (')', '}'):
                val = que.pop()
                if len(que)==0:
                    units.append(template[val:idx+1])#.strip("()\{\}"))
                    offsets.append(offset)
                    offset = 0
            elif len(que)==0:
                offset += 1

        repeat_units = []
        for unit in units:
            repeats = []
            seq = ""
            base_pattern = ''.join([c for c in unit if c not in ["(",")"]])
            for idx,i in enumerate(base_pattern):
                if i == "{":
                    repeats.append(seq)
                    seq = ""
                elif i == "}":
                    nl = []
                    for p in repeats:
                        nl.append(p+seq)
                    for p in nl:
                        repeats.append(p)
                    seq = ""
                else:
                    seq += i
            if len(seq)>0:
                repeats.append(seq)

            unwrapped = []
            for rep in repeats:
                pattern = [""]
                for char in rep:
                    if char not in tmpl.DNA_DICT:
                        for idx,p in enumerate(pattern):
                            pattern[idx] = p + char
                    else:
                        new = []
                        for iupac in tmpl.DNA_DICT[char]:
                            for p in pattern:
                                new.append(p+iupac)
                        pattern = new
                for p in pattern:
                    unwrapped.append(p)
            repeat_units.append(unwrapped)
        
        return units,repeat_units,offsets

    def get_diffs_for_all(self, template):
        """
        Gets differences between expected signals for each unit
        """
        diffs = {}
        brackets = ["(",")","{","}"]
        res = re.findall(r'[\(\{].*?[\)\}]', template)

        for r in res:
            base_pattern = ''.join([c for c in r if c not in brackets])
            pattern = [""]
            for char in base_pattern:
                if char not in tmpl.DNA_DICT:
                    for idx,p in enumerate(pattern):
                        pattern[idx] = p + char
                else:
                    new = []
                    for iupac in tmpl.DNA_DICT[char]:
                        for p in pattern:
                            new.append(p+iupac)
                    pattern = new
            for p in pattern:
                diffs[p] = self.get_consecutive_diff(p)

        return diffs

    def get_consecutive_diff(self,pattern):
        """
        Gets differences between expected signals
        """
        rep = pattern*self.kmer_size
        kmers = [rep[i:self.kmer_size+i] for i in range(len(pattern)+1)]
        pore_values = np.abs(np.diff(
            [self.pore_model[self.pore_model['kmer'] == kmer]['level_norm'].values[0] for kmer in kmers]))
        return np.mean(pore_values), np.median(pore_values)

    def collapse_repeats(self,seq):
        """
        Collapses repeat units
        """
        results = []
        slide = seq
        for i in self.repeat_units:
            if isinstance(i,list):
                results.append([0]*len(i))
            else:
                results.append(0)
        for idx,(p,off) in enumerate(zip(self.repeat_units,self.offsets)):
            slide = slide[off:]
            while len(slide)>0:
                if (isinstance(p,str) and slide[:len(p)]==p):
                    results[idx] += 1
                    slide = slide[len(p):]
                elif isinstance(p,list):
                    notfound = True
                    for idx2, possible in enumerate(p):
                        if possible==slide[:len(possible)]:
                            results[idx][idx2] += 1
                            rep = slide[len(possible):]
                            notfound = False
                else:
                    break

                if notfound:
                    break
                slide = rep
        return results

def rescale_signal(warper, state_raws, state_val, state_exp):
    """
    Rescales signal
    """
    filtered = filter_alignment(state_raws, state_val, state_exp, rescaling_config)

    filtered.sort(key = lambda x: x[0])
    filt_val = [d[0] for d in filtered]
    filt_exp = [d[1] for d in filtered]

    spline = interpolate.splrep(filt_val, filt_exp, s=len(filt_val))
    rescaled_signal = interpolate.splev(warper.template, spline)
    return rescaled_signal

def calc_ttest(arr1,arr2,win):
    """
    Calcs T-Test
    """
    std = sqrt((np.std(arr1)**2+np.std(arr2)**2)/win)
    if std==0:
        std = std + 0.0000001
    return (np.mean(arr1)-np.mean(arr2))/std


def do_segmentation(data, win):
    """
    Does segmentation
    """
    t_stats = [calc_ttest(data[idx-win:idx],data[idx:idx+win],win) for idx in range(win,len(data)-win+1)]
    segment_borders = []
    start = False
    prev = t_stats[0]

    #find peaks, i.e. segment_borders
    for idx,i in enumerate(t_stats):
        if i>3 or i<-3:
            if (i>3 and i>=prev) or (i<-3 and i<=prev):
                start = True
            else:
                if start:
                    segment_borders.append(idx+win-1)
                start = False

        elif start:
            segment_borders.append(idx+win-1)
            start = False
        prev = i

    return len(segment_borders)-1

def check_segments(segments,normal_len):
    """
    Checks segments
    """
    return [idx for idx,i in enumerate(segments) if i>=(normal_len+1)]


def dtw_wrapper(input_data):
    """
    :param input_data: list of lists with signal and reverse info
    :returns list of sequence lengths
    """

    input_signal = input_data[0]
    reverse = input_data[1]
    read_name = input_data[2]

    states = rev_states if reverse else temp_states
    endstate = rev_end_state if reverse else temp_end_state
    mask = rev_mask if reverse else temp_mask

    warper = Warper(input_signal,states,min_values_per_state,endstate,reverse)

    state_raws, state_val, state_exp = warper.create_alignment(rescaling_config['method'],rescaling_config['reps_as_one'])
    rescaled_signal = rescale_signal(warper, state_raws, state_val, state_exp)

    start, end, badmask = mask_bad_repeats(input_signal, mask, warper)
    costs = [abs(i-j) for i,j in zip(state_val[start:end],state_exp[start:end])]
    cost = np.mean(costs)

    resc_warper = Warper(rescaled_signal,states,min_values_per_state,endstate,reverse,badmask)
    resc_state_raws, resc_state_val, resc_state_exp = resc_warper.create_alignment(rescaling_config['method'],rescaling_config['reps_as_one'])

    start, end, badmask = mask_bad_repeats(input_signal, mask, resc_warper)
    resc_costs = [abs(i-j) for i,j in zip(resc_state_val[start:end],resc_state_exp[start:end])]
    resc_cost = np.mean(resc_costs)


    if visualize:
        save_warp_img(out_warp_path, FLANK_LENGTH, warper,resc_warper,reverse,read_name)

    return ((warper.seq,resc_warper.seq), (cost,resc_cost))


def mask_bad_repeats(input_signal, n_mask, warper):
    """
    Masks signal values of incorrect segments
    """
    start, end, start_idx, end_idx, bounds = find_event_borders(n_mask, warper)
    win = 3
    segment_lengths = [do_segmentation(
        input_signal[bounds[idx]-win:i+win], win) for idx, i in enumerate(bounds[1:])]
    big_segments = check_segments(segment_lengths,tr_unit)
    big_events_mask = mask_big_events(input_signal, start_idx, end_idx, bounds, big_segments)
    return start,end,big_events_mask


def find_event_borders(n_mask, warper):
    """
    Finds event borders
    """
    found_mask = [n_mask[i] for i in warper.state_transitions]
    trues = [i for i, x in enumerate(found_mask) if x]
    start,end = trues[0],trues[-1]
    start_state = warper.state_transitions[start]
    end_state = warper.state_transitions[end]
    start_idx = np.where(warper.trace==start_state)[0][0]
    end_idx = np.where(warper.trace==end_state)[0][-1]

    for_segmentation = warper.trace[start_idx:end_idx+1]

    bounds = np.where(np.diff(for_segmentation).astype(np.bool))[0]
    add = (len(bounds)-1)%tr_unit
    if add>0:
        end = end+(tr_unit-add)
        start_state = warper.state_transitions[start]
        end_state = warper.state_transitions[end]
        start_idx = np.where(warper.trace==start_state)[0][0]
        end_idx = np.where(warper.trace==end_state)[0][-1]

        for_segmentation = warper.trace[start_idx:end_idx+1]
        bounds = np.where(np.diff(for_segmentation).astype(np.bool))[0]

    bounds = start_idx + bounds
    bounds = [i for idx,i in enumerate(bounds) if idx%tr_unit==0]
    return start,end,start_idx,end_idx,bounds


def mask_big_events(input_signal, start_idx, end_idx, bounds, badones):
    """
    Masks big events
    """
    badmask = [False]*start_idx
    badmask += [False]*(bounds[0]-start_idx)

    for idx,i in enumerate(bounds[1:]):
        if idx in badones:
            badmask += [True]*(bounds[idx+1]-bounds[idx])
        else:
            badmask += [False]*(bounds[idx+1]-bounds[idx])

    badmask += [False]*(end_idx-bounds[-1])
    badmask += [False]*(len(input_signal)-end_idx)
    return badmask


def flatten_raws(arr):
    """
    Flattens list of lists
    """
    raw_indices = []
    for idx,vals in enumerate(arr):
        raw_indices.append(np.repeat(idx,len(vals)))
    raws = [val for sublist in arr for val in sublist]
    raw_indices = [val for sublist in raw_indices for val in sublist]
    return raw_indices,raws


def filter_alignment(state_raws, state_val, state_exp,config):
    """
    Filters aligned pairs of state to signal values
    """
    filtered = []
    for idx, expected_mean in enumerate(state_exp):
        if len(state_raws[idx])>=min_values_per_state and np.std(state_raws[idx])<config['max_std']:
            if abs(expected_mean - state_val[idx]) <= float(config['threshold']):
                filtered.append((state_val[idx], expected_mean))

    return filtered


def get_seqs(sequence, flanks):
    """
    Gets whole input sequence for state automaton using repeating part and flanks
    """
    tmp = flanks[0][0] + sequence + flanks[0][1]
    rev = flanks[1][0] + reverse_uniq_sequence(sequence) + flanks[1][1]
    return tmp, rev

def reverse_uniq_sequence(sequence):
    """
    Transforms the regular expression to the reverse strand
    """
    rev = ""
    for i in sequence[::-1]:
        rev += tmpl.ENCODING_DICT[i]
    return rev
