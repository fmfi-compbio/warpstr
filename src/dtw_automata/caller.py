import os
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from math import sqrt
from typing import List, Optional

import numpy as np
from scipy import interpolate

import src.dtw_automata.StateAutomata as sta
from src.config import caller_config, rescaler_config
from src.dtw_automata.plotter import save_warp_img
from src.schemas import ReadSignal
from src.squiggler.dna_sequence import get_reverse_strand


@dataclass
class StateAlignment:
    raw_values: List[float]
    expected: float
    state_value: float

    @property
    def large_enough(self) -> int:
        return len(self.raw_values) >= caller_config.min_values_per_state

    @property
    def not_outlier(self) -> bool:
        return np.std(self.raw_values) < rescaler_config.max_std

    @property
    def close(self) -> bool:
        return abs(self.expected-self.state_value) <= rescaler_config.threshold

    def good_enough(self) -> bool:
        if self.large_enough and self.not_outlier and self.close:
            return True
        else:
            return False

    @property
    def cost(self) -> float:
        return abs(self.state_value-self.expected)


@dataclass
class CallerResult:
    seq: str
    cost: float
    resc_seq: str
    resc_cost: float


@dataclass
class WarpResult:
    trace: np.ndarray

    @property
    def state_transitions(self) -> List[int]:
        return self.trace[np.insert(np.diff(self.trace).astype(np.bool8), 0, True)]

    def get_warped_signal(self, states: List[sta.State]):
        return [states[i].value for i in self.trace]

    def create_alignment(self, states: List[sta.State], signal: np.ndarray):
        alignments: List[StateAlignment] = []
        method = rescaler_config.method

        if rescaler_config.reps_as_one:
            state_indices = np.unique(self.state_transitions)
            for state_index in state_indices:
                raws = np.take(signal, np.where(self.trace == state_index)[0])
                alignments.append(
                    StateAlignment(
                        raw_values=raws,
                        state_value=process_raws(method, raws),
                        expected=states[state_index].value
                    )
                )
        else:
            idx = 0
            for state_transition in self.state_transitions:
                raws: List[float] = []
                for trace_elem in self.trace[idx:]:
                    if trace_elem != state_transition:
                        break
                    raws.append(signal[idx])
                    idx += 1
                alignments.append(
                    StateAlignment(
                        raw_values=raws,
                        state_value=process_raws(method, raws),
                        expected=states[state_transition].value
                    )
                )
        return alignments

    def delineate_tr_region(self, flank_length: int, states: List[sta.State]):
        offset = states[self.state_transitions[0]].seq_idx
        boundary = flank_length - 20
        start = bisect_left(self.trace, boundary - offset - 5)
        maxstate = self.trace[-1]
        end = bisect_right(self.trace, maxstate - boundary)
        return start, end


@dataclass
class WarpSTR:
    flank_length: int
    states: List[sta.State]
    endstate: int
    repeat_mask: List[bool]
    out_warp_path: Optional[str]
    reverse: bool
    read_name: str

    def run(self, signal: np.ndarray):

        # align STR signal with state automata
        warp_result = self.warp(signal)

        # rescale signal using the best alignment
        alignment = warp_result.create_alignment(self.states, signal)
        rescaled_signal = rescale_signal(signal, alignment)
        start, end, badmask = mask_bad_repeats(
            signal, self.repeat_mask, warp_result.trace, warp_result.state_transitions)

        # align rescaled signal with state automata
        resc_warp_result = self.warp(rescaled_signal, badmask)

        # get alignment of rescaled warp
        resc_alignment = resc_warp_result.create_alignment(self.states, rescaled_signal)
        _ = rescale_signal(rescaled_signal, resc_alignment)
        resc_start, resc_end, _ = mask_bad_repeats(
            rescaled_signal, self.repeat_mask, resc_warp_result.trace, resc_warp_result.state_transitions)

        # compute state-wise costs
        cost = np.mean([state_align.cost for state_align in alignment[start:end]])
        resc_cost = np.mean([state_align.cost for state_align in resc_alignment[resc_start:resc_end]])

        # store alignment as image for debugging
        if self.out_warp_path:
            out_path = os.path.join(
                self.out_warp_path,
                self.read_name+'_' + str(self.reverse)+'.png'
            )

            start, end = warp_result.delineate_tr_region(self.flank_length, self.states)
            resc_start, resc_end = resc_warp_result.delineate_tr_region(self.flank_length, self.states)
            if (end-start) > 2000 or (resc_end-resc_start) > 2000:
                note = f'*Truncated to 2000 as original length was {resc_end-resc_start}'
                end = start + 2000
                resc_end = resc_start + 2000
            else:
                note = None
            warped_sig = warp_result.get_warped_signal(self.states)[start:end]
            origin_sig = signal[start:end]
            resc_warped_sig = resc_warp_result.get_warped_signal(self.states)[resc_start:resc_end]
            rescaled_sig = rescaled_signal[resc_start:resc_end]
            save_warp_img(out_path, warped_sig, origin_sig, resc_warped_sig, rescaled_sig, note)

        return CallerResult(
            seq=self.get_sequence(self.reverse, warp_result),
            resc_seq=self.get_sequence(self.reverse, resc_warp_result),
            cost=cost,
            resc_cost=resc_cost
        )

    def get_sequence(self, reverse: bool, warp_result: WarpResult) -> str:
        """ Get sequence derived from traversal of states """
        seq: str = ''
        for i in warp_result.state_transitions:
            seq += self.states[i].kmer[-1]

        offset = self.states[warp_result.state_transitions[0]].seq_idx
        seq = seq[self.flank_length-offset:-self.flank_length]

        return get_reverse_strand(seq) if reverse else seq

    def warp(self, signal: np.ndarray, mask: Optional[List[bool]] = None) -> WarpResult:
        mask = mask if mask else [False]*len(signal)
        dtw_matrix = self.__calc_dtw_astates(signal, self.states, mask)
        trace = self.__backtracking(dtw_matrix, self.states, signal, mask)
        return WarpResult(trace)

    def __calc_cost(self, val1: float, val2: float):
        return abs(val1-val2)

    def __init_matrix(self, tlen: int, qlen: int):
        return np.full((tlen, qlen), np.inf)

    def __visited(self, value: np.float64):
        return False if value == np.inf else True

    def __calc_dtw_astates(self, signal: np.ndarray, states: List[sta.State], mask: List[bool]):

        dtw_matrix = self.__init_matrix(len(signal), len(states))
        # calc first state
        start_state = states[0]
        start_val = self.__calc_cost(signal[0], start_state.value)
        dtw_matrix[0, 0] = start_val

        for i in range(1, caller_config.min_values_per_state+1):
            dtw_matrix[0, i] = start_val+self.__calc_cost(signal[i], start_state.value)

        # run dtw matrix
        boundary = self.flank_length-10
        after_repeat = states[-1].seq_idx-boundary
        first_threshold = 6*boundary
        second_threshold = len(signal)-(6*boundary)
        for i, val in enumerate(signal[caller_config.min_values_per_state:], start=caller_config.min_values_per_state):
            back = caller_config.min_values_per_state - 1 if mask[i] else caller_config.min_values_per_state
            for j, curr_st in enumerate(states):
                if i < first_threshold:  # at start bound the alignment
                    if i < curr_st.seq_idx*4 and i > curr_st.seq_idx*15:
                        continue
                elif (i > second_threshold and curr_st.seq_idx < after_repeat):
                    continue

                if self.__visited(dtw_matrix[i-1, j]):
                    cost = dtw_matrix[i-1, j] + self.__calc_cost(val, curr_st.value)
                    if cost < dtw_matrix[i, j]:
                        dtw_matrix[i][j] = cost

                for prev in curr_st.incoming:
                    prev_idx = prev.idx

                    if self.__visited(dtw_matrix[i-back, prev_idx]) is False:
                        continue
                    prev_st_val = prev.value
                    skip_cost = dtw_matrix[i-back, prev_idx]

                    for skip in signal[i-back+1:i]:
                        skip_cost += self.__calc_cost(skip, prev_st_val)

                    skip_cost += self.__calc_cost(val, curr_st.value)
                    if skip_cost < dtw_matrix[i, j]:
                        dtw_matrix[i, j] = skip_cost
        return dtw_matrix

    def __backtracking(self, dtw_matrix: np.ndarray, states: List[sta.State], signal: np.ndarray, mask: List[bool]):
        """ Find the best path through the states """
        newidx = self.endstate
        last = len(dtw_matrix)-1
        trace: List[int] = []
        skip_idx = -1

        while last != 0:
            curr_dist = dtw_matrix[last][newidx]
            curr_state = states[newidx]

            if dtw_matrix[last-1][newidx] == np.inf:
                shift_delta = np.inf
            else:
                shift = dtw_matrix[last-1][newidx] + \
                    self.__calc_cost(signal[last], curr_state.value)
                shift_delta = abs(shift - curr_dist)

            skip_delta = np.inf
            if mask[last]:
                back = caller_config.min_values_per_state - 1
            else:
                back = caller_config.min_values_per_state

            for prev in curr_state.incoming:
                prev_idx = prev.idx

                # HANDLE COMING FROM NORMAL STATE
                if dtw_matrix[last-back][prev_idx] == np.inf:
                    continue
                skip_cost = dtw_matrix[last-back][prev_idx]
                for skip in signal[last-back+1:last]:
                    skip_cost += self.__calc_cost(skip, prev.value)

                skip_cost += self.__calc_cost(signal[last], curr_state.value)
                delta = abs(skip_cost-curr_dist)

                if delta < skip_delta:
                    skip_delta = delta
                    skip_idx = prev_idx

            if skip_delta < shift_delta:
                last = last - back
                trace.append(newidx)
                if skip_idx == -1:
                    raise RuntimeError('Unexpected error during backtracking')
                for _ in range(back-1):
                    trace.append(skip_idx)
                newidx = skip_idx
            else:
                last = last - 1
                trace.append(newidx)

        trace.append(newidx)
        trace.reverse()
        return np.asarray(trace, dtype=int)


def rescale_signal(signal: np.ndarray, alignments: List[StateAlignment]) -> np.ndarray:
    filtered = filter_alignment(alignments)
    filtered.sort(key=lambda x: x[0])

    filt_val = [d[0] for d in filtered]
    filt_exp = [d[1] for d in filtered]

    spline = interpolate.splrep(filt_val, filt_exp, s=len(filt_val))
    rescaled_signal = interpolate.splev(signal, spline)
    return rescaled_signal


def filter_alignment(alignments: List[StateAlignment]):
    """ Filter aligned pairs of state-to-signal values """
    return [(alignment.state_value, alignment.expected) for alignment in alignments if alignment.good_enough()]


def process_raws(method: str, vals: List[float]):
    if method == 'mean':
        return np.average(vals)
    elif method == 'median':
        return np.median(vals)
    else:
        raise KeyError('Invalid alignment method')


def mask_bad_repeats(input_signal: np.ndarray, n_mask: List[bool], trace: np.ndarray, state_transitions: List[int]):
    """
    Masks signal values of incorrect segments
    """
    start, end, start_idx, end_idx, bounds = find_event_borders(n_mask, trace, state_transitions)
    win = 3
    lengths = [segment(input_signal[bounds[idx]-win:i+win], win) for idx, i in enumerate(bounds[1:])]
    big_segments = check_segments(lengths)
    big_events_mask = mask_big_events(input_signal, start_idx, end_idx, bounds, big_segments)
    return start, end, big_events_mask


def check_segments(segment_lengths: List[int]):
    """return idxes of bigger segments """
    return [idx for idx, i in enumerate(segment_lengths) if i >= (caller_config.states_in_segment+1)]


def calc_ttest(arr1: np.ndarray, arr2: np.ndarray, win: int):
    """
    Calcs T-Test
    """
    std = sqrt((np.std(arr1)**2+np.std(arr2)**2)/win)
    if std == 0:
        std = std + 0.0000001
    return (np.mean(arr1)-np.mean(arr2))/std


def segment(data: np.ndarray, win: int):
    t_stats = [calc_ttest(data[idx-win:idx], data[idx:idx+win], win) for idx in range(win, len(data)-win+1)]
    segment_borders: List[int] = []
    start = False
    prev = t_stats[0]

    # find peaks, i.e. segment_borders
    for idx, i in enumerate(t_stats):
        if i > 3 or i < -3:
            if (i > 3 and i >= prev) or (i < -3 and i <= prev):
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


def find_event_borders(n_mask: List[bool], trace: np.ndarray, state_transitions: List[int]):
    found_mask = [n_mask[i] for i in state_transitions]
    trues = [i for i, x in enumerate(found_mask) if x]
    start, end = trues[0], trues[-1]
    start_state = state_transitions[start]
    end_state = state_transitions[end]
    start_idx = np.where(trace == start_state)[0][0]
    end_idx = np.where(trace == end_state)[0][-1]

    for_segmentation = trace[start_idx:end_idx+1]

    bounds = np.where(np.diff(for_segmentation).astype(np.bool8))[0]
    add = (len(bounds)-1) % caller_config.states_in_segment
    if add > 0:
        end = end+(caller_config.states_in_segment-add)
        start_state = state_transitions[start]
        end_state = state_transitions[end]
        start_idx = np.where(trace == start_state)[0][0]
        end_idx = np.where(trace == end_state)[0][-1]

        for_segmentation = trace[start_idx:end_idx+1]
        bounds = np.where(np.diff(for_segmentation).astype(np.bool8))[0]

    bounds = start_idx + bounds
    bounds = [i for idx, i in enumerate(bounds) if idx % caller_config.states_in_segment == 0]
    return start, end, start_idx, end_idx, bounds


def mask_big_events(input_signal: np.ndarray, start_idx: int, end_idx: int, bounds, badones):
    badmask: List[bool] = [False]*start_idx
    badmask += [False]*(bounds[0]-start_idx)

    for idx, i in enumerate(bounds[1:]):
        if idx in badones:
            badmask += [True]*(bounds[idx+1]-bounds[idx])
        else:
            badmask += [False]*(bounds[idx+1]-bounds[idx])

    badmask += [False]*(end_idx-bounds[-1])
    badmask += [False]*(len(input_signal)-end_idx)
    return badmask


# states, endstate, mask
def warpstr_call_sequential(
    input_data: ReadSignal,
    flank_length: int,
    sta: sta.StateAutomata,
    out_warp_path: Optional[str]
) -> CallerResult:
    warpstr = WarpSTR(flank_length, sta.states, sta.endstate, sta.mask,
                      out_warp_path, input_data.reverse, input_data.name)
    return warpstr.run(input_data.signal)
