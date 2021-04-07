# Copyright (C) 2021, Grzegorz Stefański - All Rights Reserved

import numpy as np
import matplotlib.pyplot as plt
from binary_graph import BinaryGraph

class NeedlemanWunch():
    '''
        NeedlemanWunch class

        Class for executing global Needleman-Wunch algorithim
        for given pair of sequences.
    '''

    def __init__(self, **kwargs):
        '''
            NeedlemanWunch class constructor

            Input:
                **kwargs: dictionary - input dictionary 

            Output:
                Object of class NeedlemanWunch
        '''
        
        self.seq1 = " " + kwargs['sequence1']
        self.seq2 = " " + kwargs['sequence2']
        self.match_score = kwargs['match_score']
        self.mismatch_score = kwargs['mismatch_score']
        self.gap_score = kwargs['gap_score']
        self.mode = kwargs['mode']
        self.print_graph = kwargs['print_graph']

    def forward(self):
        '''
            forward method

            Executing Needleman-Wunch algorithim.

            Input:
                None
            
            Output:
                None
        '''

        print(self.seq1)
        print(self.seq2)
        output = np.zeros( (len(self.seq2), len(self.seq1)) )
        s = np.zeros(3)
        path = np.array([False, False, False])

        output[:, 0] = np.arange(0, len(self.seq2) * self.gap_score, self.gap_score)
        output[0, :] = np.arange(0, len(self.seq1) * self.gap_score, self.gap_score)

        graph = BinaryGraph(self.seq1, self.seq2, retainGraph = self.print_graph)

        for i in range(1, len(self.seq2)):
            for j in range(1, len(self.seq1)):
                
                _i = i - 1
                _j = j - 1

                s[0] = output[_i, _j] + (self.match_score if self.seq2[i] == self.seq1[j] else self.mismatch_score)
                s[1] = output[i, _j] + self.gap_score
                s[2] = output[_i, j] + self.gap_score

                output[i, j] = max(s)
                path = (s == output[i, j])

                if (i == 1 or j == 1) and path[0]:
                    graph.addRoot(i, j)
                
                elif path[0]:
                    graph.addNode(i, j, _i, _j, self.seq1[j], self.seq2[i])

                if path[1]:
                    graph.addNode(i, j, i, _j, self.seq1[j], "-")
                
                if path[2]:
                    graph.addNode(i, j, _i, j, "-", self.seq2[i])


        graph.scores = output
        print()
        print(output)
        # self.plot(output)

        print("\nRoots")
        print(graph.printRoots())

        print("\nNodes")
        print(graph.printNodes())

        
        if self.print_graph:
            print()
            # graph.printTree()

        print()
        graph.getPaths(self.mode, self.gap_score)

    def plot(self, matrix):
        '''
            plot methods

            Method for ploting computed score matrix

            Input:
                matrix: 2d vector - score matrix

            Output:
                None
        '''


        labels_x = list(self.seq1)
        labels_y = list(self.seq2)

        fig, ax = plt.subplots()
        im = ax.imshow(matrix)

        # We want to show all ticks...
        ax.set_xticks(np.arange(len(labels_x)))
        ax.set_yticks(np.arange(len(labels_y)))
        # ... and label them with the respective list entries
        ax.set_xticklabels(labels_x)
        ax.set_yticklabels(labels_y)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(labels_y)):
            for j in range(len(labels_x)):
                text = ax.text(j, i, matrix[i, j],
                            ha="center", va="center", color="w")

        ax.set_title("Harvest of local farmers (in tons/year)")
        fig.tight_layout()
        plt.show()

        plt.waitforbuttonpress()