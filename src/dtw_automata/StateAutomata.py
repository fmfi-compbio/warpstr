import src.templates as tmpl


class StateAutomata:
    """
    Class State Automata
    """

    def __init__(self, pore_model, sequence, state_size):
        self.verbose = 0
        self.sequence = sequence
        self.state_size = state_size
        self.pore_model = pore_model
        self.endstate = -1
        self.repstart = -1
        self.repend = -1
        self.simple_states = self.create_simple_states()
        self.kmer_states, self.mask = self.create_kmer_states()

    def set_repstart(self, idx):
        """
        Sets the start of repeating part
        """
        if self.repstart == -1:
            self.repstart = idx

    def set_repend(self, idx):
        """
        Sets the end of repeating part
        """
        self.repend = idx

    def create_simple_states(self):
        """
        Creates simple state automata from input sequence
        :param input_seq: formatted sequence
        :return simple_states: list of created simple states
        """

        simple_states = []

        #serves purpose for handling repeats, que for () and wque for {}
        que, wque = [], []
        spec_states, prev_states = [], []

        #initialize first simple state
        prev_states.append(Simple_state(self.sequence[0], 0))
        simple_states.append(prev_states[0])
        idx = 1
        for i in self.sequence[1:]:

            #store this index in queue so later it can be added when ) is found
            if i == '(':
                que.append(idx)
                self.set_repstart(idx)
            #add path from last simple state to the simple state with index from que
            #which represents when this repeat started with (
            elif i == ')':
                self.set_repend(idx)
                last_idx = que.pop()

                for prev in prev_states:
                    if isinstance(last_idx, list):
                        for last_idx_iter in last_idx:
                            prev.addnextpos(simple_states[last_idx_iter])
                    else:
                        prev.addnextpos(simple_states[last_idx])

            #similarly as for ( but there can be zero loops
            elif i == '{':
                wque.append(idx-1)
                #spec_states.append(simple_states[last_idx])
            elif i == '}':
                last_idx = wque.pop()
                #for p in prev_states:
                #    if type(last_idx) is list:
                #        for last_idx_iter in last_idx:
                #            p.addnextpos(simple_states[last_idx_iter+1])
                #    else:
                #        p.addnextpos(simple_states[last_idx+1])
                spec_states.append(simple_states[last_idx])

            #decode additional DNA alphabet characters
            elif i in tmpl.DNA_DICT.keys():

                #creates list of new states
                new_states = []
                if (len(que) > 0 and idx == que[-1]) or (len(wque) > 0 and idx == wque[-1]):
                    check = True
                else:
                    check = False

                nque = []
                for j in tmpl.DNA_DICT[i]:
                    new_state = Simple_state(j, idx)
                    simple_states.append(new_state)

                    if check:
                        nque.append(idx)
                    idx += 1
                    new_states.append(new_state)

                    #add chain from previous states to the new state
                    for prev in prev_states:
                        prev.addnextpos(new_state)

                    #for p in spec_states: #SHOULD I?
                    #    p.addnextpos(new_state) #SHOULD I?
                if check:
                    del que[-1]
                    que.append(nque)
                #list of previous states is updated with the list of new states
                prev_states = new_states
                #spec_states = [] #SHOULD I?

            #handle basic DNA character
            else:
                #print("loaded normal")
                new_state = Simple_state(i, idx)
                simple_states.append(new_state)

                #add chain from previous states to the new state
                for prev in prev_states:
                    prev.addnextpos(new_state)
                prev_states = []

                #list of previous states is updated with the new state
                prev_states.append(new_state)
                for prev in spec_states:
                    prev.addnextpos(new_state)
                spec_states = []
                idx += 1

        return simple_states

    ###KMER STATES
    def create_kmer_states(self):
        """
        Creates k-mer state automata from simple states automata
        :param simple_states: list of simple states
        :param pore_model: loaded normalized pore model table
        :param STATE_SIZE: size of k-mer
        """

        #initialize empty list for each possible simple state
        kmer_states = []
        for i in range(len(self.simple_states)):
            kmer_states.append([])

        #prepare first k-mer state
        curr_kmer = ""
        idx = 0
        while len(curr_kmer) != self.state_size:
            curr_kmer += self.simple_states[idx].kmer
            idx = idx + 1
        idx = idx - 1
        curr_state = self.simple_states[idx]
        new_kstate = Kmer_state(curr_kmer, len(
            kmer_states[idx]), idx, curr_state.nextpos, -1, -1)
        kmer_states[idx].append(new_kstate)

        #new k-mer state is added to the queue for processing
        que = []
        que.append(new_kstate)

        #while queue is not empty, pop k-mer state and create next states
        while len(que) > 0:
            curr_kmer_state = que[-1]
            del que[-1]

            curr_kmer = curr_kmer_state.kmer[1:]
            for next_simp_st in curr_kmer_state.to_visit:
                new_kmer = curr_kmer + next_simp_st.kmer
                idx = next_simp_st.idx
                dontstop = True
                for test in kmer_states[idx]:
                    if new_kmer == test.kmer and checkprevious(curr_kmer_state, test):
                        dontstop = False
                        curr_kmer_state.chain(test)
                if dontstop:
                    newidx = len(kmer_states[idx])
                    new_kstate = Kmer_state(new_kmer, newidx, idx, next_simp_st.nextpos,
                                            curr_kmer_state.simple_idx, curr_kmer_state.idx)
                    curr_kmer_state.chain(new_kstate)
                    kmer_states[idx].append(new_kstate)
                    que.append(new_kstate)

        #flatten list of lists of kmer_states to the basic 1D list
        #while also setting k-mer value from pore model table
        repeat_mask = []
        states = []
        for i in kmer_states:
            for j in i:
                j.set_realidx(len(states))
                states.append(j)
                if j.simple_idx >= self.repstart-1 and j.simple_idx <= self.repend+10:
                    repeat_mask.append(True)
                else:
                    repeat_mask.append(False)
                j.set_value(self.pore_model)

        for j in i:
            self.endstate = j.realidx

        if self.verbose == 1:
            prettyprint(kmer_states)

        for i in states:
            for j in i.nextpos:
                j.incoming.append(i)

        return states, repeat_mask


class Simple_state:
    """
    Class for simple one nucleotide states.
    Each nucleotide of input sequence is stored separately here.
    Attribute nextpos tells which simple states can follow from this state
    """

    def __init__(self, kmer, idx):
        """
        :param kmer: nucleotide representing state
        :param idx: tells position of the input nucleotide in sequence
        """
        self.kmer = kmer
        self.nextpos = []
        self.idx = idx

    def addnextpos(self, nextstate):
        """
        Adds state position as next position of that state
        """
        self.nextpos.append(nextstate)


class Kmer_state:
    """
    Class for k-mer states.
    Using simple states possible k-mer states are derived and stored in this state
    Attribute nextpos tells which simple states can follow from this state
    """

    def __init__(self, kmer, idx, simple_idx, nextstates, prevsimple, previdx):
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
        self.nextpos = []
        self.idx = idx
        self.simple_idx = simple_idx
        self.to_visit = nextstates
        self.previous = prevsimple, previdx
        self.realidx = -1
        self.value = -1
        self.incoming = []

    def chain(self, state):
        """
        :param state: chains current k-mer state with another by adding into the nextpos list
        """
        self.nextpos.append(state)

    def print_state(self):
        """
        Prints state information
        """
        next_simple_idx = [i.idx for i in self.to_visit]
        next_idx = [(i.simple_idx, i.idx) for i in self.nextpos]
        print("State", "{:03d}".format(self.simple_idx), "kmer:", self.kmer,
              "visit:", next_simple_idx, "next:", next_idx, "prev:", self.previous)
        #print("")

    def set_realidx(self, realidx):
        """
        Setter of realidx
        """
        self.realidx = realidx

    def set_value(self, kmertable):
        """
        Setter of state value from kmer table
        """
        self.value = kmertable[kmertable['kmer']
                               == self.kmer]['level_norm'].values[0]


def checkprevious(curr_state, next_state):
    """
    Check whether a loop of the same states is not formed
    """
    to_be_added = (curr_state.simple_idx, curr_state.idx)
    existing = next_state.previous
    if to_be_added[0] == existing[0]:
        return True
    else:
        return False
