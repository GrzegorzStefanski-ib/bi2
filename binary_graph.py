# Copyright (C) 2021, Grzegorz Stefański - All Rights Reserved

import numpy as np

class BinaryGraph():
    '''
        BinaryGraph class

        Data structure class. Specific type of binary tree 
        where every node can have multiple parents as well 
        as childrens.
    '''

    def __init__(self, seq1, seq2, retainGraph):
        '''
            Constructor of BinaryGraph class

            Input: 
                seq1: string - Sequence 1 string.
                seq2: string - Sequence 2 string.

            Output:
                BinaryGraph: Constructed object of class BinaryGraph
        '''

        self.seq1 = seq1
        self.seq2 = seq2

        self.scores = []

        self.roots = {}
        self.leafs = []
        self.nodes = {}

        self.retainGraph = retainGraph

    def addLeaf(self, x, y):
        '''
            addLeaf method

            Gets node with (x, y) coordinates and adds it to the 
            list of leafs (nodes without parents).  

            Input:
                x: int - vertical coordinate of sequences matrix
                y: int - horizontal coordinate of sequences matrix

            Output:
                None
        '''

        coordinates = str(x) + "," + str(y)
        if coordinates in self.nodes:
            self.leafs.append(self.nodes[coordinates])

    def addRoot(self, x, y):
        '''
            addRoot method

            Adds new root node and creates padding nodes if required.

            Input:
                x: int - vertical coordinate of sequences matrix
                y: int - horizontal coordinate of sequences matrix

            Output:
                None
        '''

        node = None
        s1 = ''
        s2 = ''

        if x != 1:
            if self.retainGraph:
                node = Node(self, self, 1, 0)
                node.addLetters("-", self.seq2[1])

                if(node.getCoordinates not in self.roots):
                    self.roots[node.getCoordinates()] = node

                else:
                    node = self.roots[node.getCoordinates]

                for it in range(2, x):
                    node = self.addNode(it, 0, it - 1, 0, "-", self.seq2[it])
            
            else:
                for it in range(1, x):
                    s1 += "-"
                    s2 += self.seq2[it]
        
        elif y != 1:
            if self.retainGraph:
                node = Node(self, self, 0, 1)
                node.addLetters(self.seq1[1], "-")

                if(node.getCoordinates not in self.roots):
                    self.roots[node.getCoordinates()] = node

                else:
                    node = self.roots[node.getCoordinates]

                for it in range(2, y):
                    node = self.addNode(0, it, 0, it - 1, self.seq1[it], "-")

            else:
                for it in range(1, y):
                    s1 += self.seq1[it]
                    s2 += "-"

        if node == None:        
            node = Node(self, self, x, y)
            node.addLetters(s1 + self.seq1[y], s2 + self.seq2[x])

            if(node.getCoordinates not in self.roots):
                self.roots[node.getCoordinates()] = node

            else:
                node = self.roots[node.getCoordinates]

        else: 
            node = self.addNode(x, y, node.x, node.y, self.seq1[y], self.seq2[x])

    def registerNode(self, node):
        '''
            registerNode method

            Registers node in graph if does not exist.

            Input:
                node - node to register.

            Output:
                None
        '''
        
        coordinates = node.getCoordinates()
        
        if coordinates not in self.nodes:
            self.nodes[coordinates] = node

    
    def addNode(self, node_x, node_y, parent_x, parent_y, letter_x, letter_y):
        '''
            addNode method

            Adds new node to graph.

            Input:
                node_x: int - node vertical coordinate in sequences matrix.
                node_y: int - node horizontal coordinate in sequences matrix.
                parent_x: int - parent vertical coordinate in sequences matrix.
                parent_y: int - parent horizontal coordinate in sequences matrix.
                letter_x: char - letter corresponding to node position in sequence1.
                letter_y: char - letter corresponding to node position in sequence2.

            Output:
                Node: Added node
        '''

        parents = []
        
        parent_coordinates = str(parent_x) + "," + str(parent_y)
        
        if parent_coordinates in self.roots:
            parents.append(self.roots[parent_coordinates])

        if parent_coordinates in self.nodes:
            parents.append(self.nodes[parent_coordinates])

        node_coordinates = str(node_x) + "," + str(node_y)

        if node_coordinates in self.nodes:
            node = self.nodes[node_coordinates]
        else:
            node = Node(parents[0], self, node_x, node_y)

        for parent in parents:
            for it in range( len(parent.letter_x) ):
                node.addLetters(parent.letter_x[it] + letter_x, parent.letter_y[it] + letter_y)
            
            parent.regiserChild(node)

        if not self.retainGraph:
            diagonal_coordinates = str(node_x - 1) + "," + str(node_y - 1)
            
            if diagonal_coordinates in self.nodes and len(self.nodes[diagonal_coordinates].childrens) != 0:
                del self.nodes[diagonal_coordinates]

        return node

    def getPaths(self, mode, gap_score):
        '''
            getPaths method

            Finds and prints determinated paths in respect to the specified mode.

            Input:
                mode: string - Mode for path finding. One of:
                    "all" - prints all paths.
                    "full_paths" - prints only paths that covers both sequences fully. 
                    "top_score" - prints only full paths with best scores.
                gap_score: float - score for gap.

            Output:
                None
        '''

        l1 = len( self.seq1 ) - 1
        l2 = len( self.seq2 ) - 1
        
        s = ""

        if mode == "top_score":
            
            max_score = None
            max_coordinates = []

            for it in range(0, l1 + 1):
                sc = self.scores[l2, it] + gap_score * (l1 - it)
                
                if max_score == None or sc > max_score:
                    max_score = sc
                    max_coordinates = [[l2, it]]

                elif sc == max_score:
                    max_coordinates += [[l2, it]]

            for it in range(0, l2):
                sc = self.scores[it, l1] + gap_score * (l2 - it)

                if max_score == None or sc > max_score:
                    max_score = sc
                    max_coordinates = [[it, l1]]

                elif sc == max_score:
                    max_coordinates += [[it, l1]]                

        for node in self.nodes.values():
            if len(node.childrens) == 0:
                for it in range( len(node.letter_x) ):
                    
                    if  ( mode == "full_path"  and (l1 == len(node.letter_x[it].replace("-", "")) or l2 == len(node.letter_y[it].replace("-", "")) )) or \
                        ( mode == "top_score" and [node.x, node.y] in max_coordinates and (l1 == len(node.letter_x[it].replace("-", "")) or l2 == len(node.letter_y[it].replace("-", "")) )) or \
                        ( mode == "all" ):
                            s += node.letter_x[it] + "\n" 
                            s += node.letter_y[it] + "\n"
                            s += "score: " + str(self.scores[node.x, node.y]) + "\n\n"
                        
        print(s)

    def printTree(self):
        '''
            printTree method

            Prints graph structure in form of binary tree (for testing purpose).
            If graph is not retained then it will only print last nodes as 
            information about other nodes is not being stored.

            Input:
                None

            Output:
                None
        '''

        print("BinaryTree")

        l = len(self.roots) - 1
        i = 0

        for root in self.roots.values():
            root.printNode(" └" if i == l else " ├")
            i += 1


    def printRoots(self):
        '''
            printRoots method

            Prints roots cooridnates and letters (for testing purpose).

            Input:
                None

            Output:
                None
        '''

        for root in self.roots.values():
            print(root.x, root.y, root.letter_x, root.letter_y)
    
    def printNodes(self):
        '''
            printNodes method

            Prints nodes cooridnates and letters (for testing purpose).

            Input:
                None

            Output:
                None
        '''

        for nodes in self.nodes.values():
            print(nodes.x, nodes.y, nodes.letter_x, nodes.letter_y)

class Node():

    '''
        Node class

        Node for binary graph. It stores information about
        node coordinates (in sequences plane), corresponding 
        letters, and childrens.
    '''

    def __init__(self, parent, tree, x, y):
        '''
            Constructor of Node class

            Input: 
                parent: Node - parent node.
                tree: BinaryGraph - pointer to main data structure.
                x: int - horizontal coordinate of node
                y: int - vertical coordinate of node

            Output:
                Node: Constructed object of class Node
        '''
        
        self.tree = tree
        self.parent = parent

        self.x = x
        self.y = y

        self.letter_x = []
        self.letter_y = []

        self.childrens = []

    def addLetters(self, letter_x, letter_y):
        '''
            Adds path to node

            Input: 
                letter_x: string - path to add for sequence 1
                letter_y: string - path to add for sequence 2

            Output:
                None
        '''

        if (letter_x not in self.letter_x) or (letter_y not in self.letter_y):
            self.letter_x.append(letter_x)
            self.letter_y.append(letter_y)

    def regiserChild(self, node):
        '''
            Registers new child of node.

            Input: 
                node: Node - child node to register

            Output:
                None
        '''

        self.childrens.append( node )
        self.tree.registerNode( node )

    def printNode(self, heading):
        '''
            Prints node informations to console.

            Input: 
                heading: string - heading (prefix) to add to printed string 

            Output:
                None
        '''

        print(heading, "Node", self.x, self.y, self.letter_x, self.letter_y)

        l = len(self.childrens)

        heading = heading.replace("└", "  ").replace("├", "│ ")

        for it in range(l):
            if(it != l - 1):
                self.childrens[it].printNode(heading + "├")
            else:
                self.childrens[it].printNode(heading + "└")

    def getCoordinates(self):
        '''
            Returns node coordinates string.

            Input: 
                None

            Output:
                string: coordinates string (used in BinaryGraph for indexing node dictionary)
        '''

        return str(self.x) + "," + str(self.y)