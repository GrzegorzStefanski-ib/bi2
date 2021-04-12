# Copyright (C) 2021, Grzegorz Stefański - All Rights Reserved

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


from binary_graph import BinaryGraph

class NeedlemanWunch():
    '''
        NeedlemanWunch class

        Class for executing global Needleman-Wunch algorithm
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

            Executing Needleman-Wunch algorithm.

            Input:
                None
            
            Output:
                None
        '''

        print()
        print("Sequence 1:", self.seq1)
        print("Sequence 2:", self.seq2)
        print()

        output = np.zeros( (len(self.seq2), len(self.seq1)) )
        s = np.zeros(3)

        output[:, 0] = np.arange(0, len(self.seq2) * self.gap_score, self.gap_score)
        output[0, :] = np.arange(0, len(self.seq1) * self.gap_score, self.gap_score)

        graph = BinaryGraph(self.seq1, self.seq2, retainGraph = self.print_graph)

        for i in tqdm(range(1, len(self.seq2))):
            for j in range(1, len(self.seq1)):
                
                _i = i - 1
                _j = j - 1

                s[0] = output[_i, _j] + (self.match_score if self.seq2[i] == self.seq1[j] else self.mismatch_score)
                s[1] = output[i, _j] + self.gap_score
                s[2] = output[_i, j] + self.gap_score

                output[i, j] = max(s)
                path = (s == output[i, j])

                if output [i,j]> min(i,j) * -3:

                    if (i == 1 or j == 1) and path[0]:
                        graph.addRoot(i, j)
                        
                    elif path[0]:
                        graph.addNode(i, j, _i, _j, self.seq1[j], self.seq2[i])

                    if path[1]:
                        graph.addNode(i, j, i, _j, self.seq1[j], "-")
                        
                    if path[2]:
                        graph.addNode(i, j, _i, j, "-", self.seq2[i])



        graph.scores = output

        if self.print_graph:
            print()
            graph.printTree()

        print()
        path = graph.getPaths(self.mode, self.gap_score)


        self.plot(output, path)
        self.save_to_file(path)

    def plot(self, matrix, path):
        '''
            plot methods

            Method for ploting computed score matrix

            Input:
                matrix: 2d vector - score matrix
                path: list (2, n) - optimal path

            Output:
                None
        '''

        x = [0]
        y = [0]
        
        for it in range( len( path[0] )):
            
            if path[0][it] != "-":
                x.append( x[-1] + 1)
                x[-2] += 0.25
            else:
                x.append( x[-1] )

            if path[1][it] != "-":
                y.append( y[-1] + 1 )
                y[-2] += 0.25
            else:
                y.append( y[-1] )

        labels_x = list(self.seq1)
        labels_y = list(self.seq2)

        fig, ax = plt.subplots(figsize=(min(16, 6 + len(self.seq1) // 15) , min(16, 6 + len(self.seq2) // 15) ))
        im = ax.imshow(matrix)


        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 

        if len(labels_x) < 30 and len(labels_y) < 30:

            ax.set_xticks( np.arange( len(labels_x) ))
            ax.set_yticks( np.arange( len(labels_y) ))

            ax.set_xticklabels(labels_x)
            ax.set_yticklabels(labels_y)

            ax.set_xticklabels(labels_x)
            ax.set_yticklabels(labels_y)

            for i in range(len(labels_y)):
                for j in range(len(labels_x)):
                    text = ax.text(j, i, matrix[i, j],
                                ha="center", va="center", color="w")

        else:
            ax.set_xlabel(self.seq1)
            ax.set_ylabel(self.seq2)

        plt.scatter(x, y, c = "black", s = 20 / (len( path[0] ) // 15 + 1))

        fig.tight_layout()
        plt.show()

    def save_to_file(self, path):

        l = len( path[0] )

        score = 0
        match = 0
        gaps = 0
        mismach = 0

        for it in range(l):
            
            if path[0][it] != "-" and path[1][it] != "-":
                
                if path[0][it] !=  path[1][it]:
                    score += self.mismatch_score
                    mismach += 1
                
                else:
                    score += self.match_score
                    match += 1

            else:
                score += self.gap_score
                gaps += 1

        s = "Sequence 1:" + self.seq1 + "\nSequence 2:" + self.seq2 + "\nMatch: " + str(match)
        s += "\nMismatch: " + str(mismach) + "\nGap: " + str(gaps)
        s += "\nScore: " + str(score) + "\nLength: " + str(l)
        s += "\nIdentity: " + str(match) + "/" + str(l) + " (" + str(match * 100 // l) + "%)"
        s += "\nGaps: " + str(gaps) + "/" + str(l) + " (" + str(gaps * 100 // l) + "%)"

        f = open("out.txt", "w")
        f.write(s)
        f.close()

