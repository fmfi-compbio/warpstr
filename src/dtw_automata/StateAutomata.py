from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Union

import src.templates as tmpl
from src.squiggler.pore_model import pore_model


@dataclass
class State:
    kmer: str
    value: float
    seq_idx: int
    idx: int
    incoming: List[State] = field(default_factory=list)
    nextpos: List[int] = field(default_factory=list)


@dataclass
class SimpleState:
    """
    Class for simple one nucleotide states.
    Each nucleotide of input sequence is stored separately here.
    Attribute nextpos tells which simple states can follow from this state
    """
    kmer: str
    idx: int
    nextpos: List[int] = field(default_factory=list)

    def addnextpos(self, nextstate_idx: int):
        self.nextpos.append(nextstate_idx)


@dataclass
class StateAutomata:
    states: List[State]
    endstate: int
    mask: List[bool]
    repstart: int = -1
    repend: int = -1

    def __init__(self, sequence: str):
        simple_states = self._create_simple_states(sequence)
        states, mask, end = self._create_kmer_states(simple_states)
        self.states = states
        self.mask = mask
        self.endstate = end

    def _set_repstart(self, idx: int):
        if self.repstart == -1:
            self.repstart = idx

    def _set_repend(self, idx: int):
        self.repend = idx

    def _create_simple_states(self, sequence: str):
        """
        Creates simple state automata from input sequence
        :param input_seq: formatted sequence
        :return simple_states: list of created simple states
        """

        simple_states: List[SimpleState] = []

        # serves purpose for handling repeats, que for () and wque for {}
        que: Union[List[int], List[List[int]]] = []
        wque: Union[List[int], List[List[int]]] = []
        spec_states: List[SimpleState] = []
        prev_states: List[SimpleState] = []

        # initialize first simple state
        root = SimpleState(sequence[0], 0)
        prev_states.append(root)
        simple_states.append(root)
        idx = 1
        for i in sequence[1:]:

            # store this index in queue so later it can be added when ) is found
            if i == '(':
                que.append(idx)
                self._set_repstart(idx)
            # add path from last simple state to the simple state with index from que
            # which represents when this repeat started with (
            elif i == ')':
                self._set_repend(idx)
                last_idx = que.pop()

                for prev in prev_states:
                    if isinstance(last_idx, list):
                        for last_idx_iter in last_idx:
                            prev.addnextpos(simple_states[last_idx_iter].idx)
                    else:
                        prev.addnextpos(simple_states[last_idx].idx)

            # similarly as for ( but there can be zero loops
            elif i == '{':
                wque.append(idx-1)
            elif i == '}':
                last_idx = wque.pop()
                spec_states.append(simple_states[last_idx])

            # decode additional DNA alphabet characters
            elif i in tmpl.DNA_DICT.keys():

                # creates list of new states
                new_states: List[SimpleState] = []
                if (len(que) > 0 and idx == que[-1]) or (len(wque) > 0 and idx == wque[-1]):
                    check = True
                else:
                    check = False

                nque: List[int] = []
                for j in tmpl.DNA_DICT[i]:
                    new_state = SimpleState(j, idx)
                    simple_states.append(new_state)

                    if check:
                        nque.append(idx)
                    idx += 1
                    new_states.append(new_state)

                    # add chain from previous states to the new state
                    for prev in prev_states:
                        prev.addnextpos(new_state.idx)

                if check:
                    del que[-1]
                    que.append(nque)
                # list of previous states is updated with the list of new states
                prev_states = new_states

            # handle basic DNA character
            else:
                new_state = SimpleState(i, idx)
                simple_states.append(new_state)

                # add chain from previous states to the new state
                for prev in prev_states:
                    prev.addnextpos(new_state.idx)
                prev_states = []

                # list of previous states is updated with the new state
                prev_states.append(new_state)
                for prev in spec_states:
                    prev.addnextpos(new_state.idx)
                spec_states = []
                idx += 1

        return simple_states

    def _create_kmer_states(self, simple_states: List[SimpleState]):

        # initialize empty list for each possible simple state
        kmer_states: List[List[KmerState]] = []
        for i in range(len(simple_states)):
            kmer_states.append([])

        # prepare first k-mer state
        curr_kmer = ''
        idx = 0
        while len(curr_kmer) != pore_model.kmersize:
            curr_kmer += simple_states[idx].kmer
            idx = idx + 1
        idx = idx - 1
        curr_state = simple_states[idx]
        new_kstate = KmerState(curr_kmer, len(kmer_states[idx]), idx, curr_state.nextpos, -1, -1)
        kmer_states[idx].append(new_kstate)

        # new k-mer state is added to the queue for processing
        que: List[KmerState] = []
        que.append(new_kstate)

        # while queue is not empty, pop k-mer state and create next states
        while len(que) > 0:
            curr_kmer_state = que[-1]
            del que[-1]

            curr_kmer = curr_kmer_state.kmer[1:]
            for next_simp_st_idx in curr_kmer_state.to_visit_idx:
                next_simp_st = simple_states[next_simp_st_idx]
                new_kmer = curr_kmer + next_simp_st.kmer
                idx = next_simp_st.idx
                dontstop = True
                for test in kmer_states[idx]:
                    if new_kmer == test.kmer and checkprevious(curr_kmer_state, test):
                        dontstop = False
                        curr_kmer_state.chain(test)
                if dontstop:
                    newidx = len(kmer_states[idx])
                    new_kstate = KmerState(new_kmer, newidx, idx, next_simp_st.nextpos,
                                           curr_kmer_state.simple_idx, curr_kmer_state.idx)
                    curr_kmer_state.chain(new_kstate)
                    kmer_states[idx].append(new_kstate)
                    que.append(new_kstate)

        # flatten list of lists of kmer_states to the basic 1D list
        # while also setting k-mer value from pore model table
        repeat_mask: List[bool] = []
        counter = 0
        for i in kmer_states:
            for j in i:
                j.realidx = counter
                counter = counter + 1
                # states.append(j.simplify(len(states)))
                if j.simple_idx >= self.repstart-1 and j.simple_idx <= self.repend+10:
                    repeat_mask.append(True)
                else:
                    repeat_mask.append(False)

        i = kmer_states[-1]
        endstate = -1
        for j in i:
            endstate = j.realidx

        states: List[State] = []
        for i in kmer_states:
            for j in i:
                states.append(j.simplify(len(states)))

        for i in states:
            for j in i.nextpos:
                next_st = states[j]
                next_st.incoming.append(i)

        return states, repeat_mask, endstate


class KmerState:
    """
    Class for k-mer states.
    Using simple states possible k-mer states are derived and stored in this state
    Attribute nextpos tells which simple states can follow from this state
    """
    realidx: int

    def __init__(self, kmer: str, idx: int, simple_idx: int, nextstates: List[int], prevsimple: int, previdx: int):
        """
        :param kmer: k-mer representing state
        :param nexpost: possible next k-mer states
        :param idx: tells the number of k-mer states for simple_idx
        :param simple_idx: tells position of the input nucleotide in sequence
        :param to_visit: list of states to visit next
        :param previous: tells from which index of simple state and k-mer state we got here
        :param realidx: stores true index of this object in list of k-mer states
        :param value: value of k-mer state from the pore model table
        """
        self.kmer = kmer
        self.nextpos: List[KmerState] = []
        self.idx = idx
        self.simple_idx = simple_idx
        self.to_visit_idx = nextstates
        self.previous = prevsimple, previdx

    def simplify(self, idx: int):
        return State(
            kmer=self.kmer,
            nextpos=[kmerstate.realidx for kmerstate in self.nextpos],
            value=pore_model.get_value(self.kmer),
            seq_idx=self.simple_idx,
            idx=idx
        )

    def chain(self, state: KmerState):
        """
        :param state: chains current k-mer state with another by adding into the nextpos list
        """
        self.nextpos.append(state)


def checkprevious(curr_state: KmerState, next_state: KmerState) -> bool:
    """
    Check whether a loop of the same states is not formed
    """
    to_be_added = (curr_state.simple_idx, curr_state.idx)
    existing = next_state.previous
    if to_be_added[0] == existing[0]:
        return True
    else:
        return False
