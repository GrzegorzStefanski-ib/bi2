import sys
from os import path
import re
import click

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from needleman_wunch import NeedlemanWunch
import NCBI

sequence_regex = "^[ACDEFGHIKLMNPQRSTVWY\s]+$"


def normalize_sequence(sequence):
    '''
        normalize_sequence

        Normalizes sequence string to expected format.

        input:
            sequence: string - string with sequence to normalize.
        
        output:
            string: Normalized sequence
    '''

    return sequence.upper().replace(" ", "").replace("\n", "")

def prepare_sequence(sequence):
    '''
        prepare_sequence

        Preperes sequence for computing. If input is direct sequence it is normalized.
        If input is path to .fasta file it attempt to read it and then normalize included sequence.
        If input is id it attempts to download .fasta file from NCBI and then read it.

        input:
            sequence: string - string with sequence/path to .fasta file/id for NCBI to read.
        
        output:
            string: Sequence
    '''
    if re.search(sequence_regex, sequence, re.IGNORECASE):
        return normalize_sequence(sequence)

    elif path.isfile(sequence) or path.isfile(sequence + ".fasta") or path.isfile(sequence + ".FASTA"):
        return normalize_sequence( fasta_read(sequence)[0] )

    else:
        filepath = NCBI.download(sequence)

        if filepath == None or not path.isfile(filepath):
            raise click.BadParameter('')

        return normalize_sequence( fasta_read(filepath)[0] )        


def fasta_read(directory):
    '''
        fasta_read

        Reads .fasta file and returns vector of sequences from file.

        input:
            directory: string - path to .fasta file.

        output:
            string [n]: Vector of n sequences from file.
    '''

    if not ( directory.endswith(".fasta") or directory.endswith(".FASTA") ):
        directory += ".fasta"

    sequences = []
    seq = ""

    with open(directory, "r") as file_handle:
        for line in file_handle.readlines():
            
            if line[0] == ">":
                if len(seq) != 0:
                    sequences.append(seq.upper())       
                seq = ""
        
            else:
                seq += line
    
    if len(seq) != 0:
        sequences.append(seq.upper())

    if len(sequences) == 0:
        raise click.BadParameter("File: " + directory + " does not contain any sequence or is not in right format (.fasta).")
        
    return sequences

def validate_sequence(ctx, param, value):
    '''
        validate_sequence

        Validates sequence/.fasta file/NCBI id and prepare sequence for computing 
        
        input:
            ctx: click.core.Context - context object.
            param: string - parameter name.
            value: string - sequence.

        output:
            string: Sequence.
    '''
    
    try:
        return prepare_sequence(value)

    except ValueError:
        raise click.BadParameter('')

def validate_score(ctx, param, value):
    '''
        validate_score

        Validates score  
        
        input:
            ctx: click.core.Context - context object.
            param: string - parameter name.
            value: int - score.

        output:
            int: Score.
    '''

    if not isinstance(value, int):
        raise click.BadParameter('')
    
    return value

@click.command()

@click.option('--sequence1', '--seq1', '--s1', prompt = 'Sequence 1', callback = validate_sequence, help = 'Sequence 1')
@click.option('--sequence2', '--seq2', '--s2', prompt = 'Sequence 2', callback = validate_sequence, help = 'Sequence 2')

@click.option('--match_score', '--ms',     prompt = 'Match score',    help = 'Match score',    type = int, callback = validate_score)
@click.option('--mismatch_score', '--mms', prompt = 'Mismatch score', help = 'Mismatch score', type = int, callback = validate_score)
@click.option('--gap_score', '--gs',       prompt = 'Gap score',      help = 'Gap score',      type = int, callback = validate_score)

@click.option('--mode', '--m', default = "top_score", prompt = 'Mode', show_default = True, type = click.Choice(['all', 'full_path', 'top_score'], case_sensitive=False), help = 'Result filtering mode.')
@click.option('--print_graph', '--pg', default = False, prompt = 'Print graph', show_default = True, type=bool,  help = 'Print constructed graph.')

def main(**kwargs):

    '''
        Algorytm Needleman'a-Wunch'a 

        Examples:

         \b
         # Longer version
         python main.py --sequence1=AAA --sequence2=CCC --match_score=1 --mismatch_score=-1 --gap_score=-2 --mode=top_score --print_graph=False

         \b
         # Shorter version
         python main.py --s1=AAA --s2=CCC --ms=1 --mms=-1 --gs=-2 --m=top_score --pg=False

         \b
         # You can also run this script without any arguments.
         python main.py

         \b
         Mode (--mode):
          all          Display all possible paths.
          full_path    Display only results that covers whole sequences.
          top_score    Display only paths with higher score.

    '''

    # print(kwargs)
    # ali = Aligator()

    # Aligator("TTACG", "TTCG", 1, -0.5, -1).forward()
    # Aligator("TTCG", "TTACG", 1, -0.5, -1).forward()
    # Aligator("TTACG", "TTCG", 1, -1, -2).forward()
    # Aligator("GAAC", "CAAGAC", 1, -1, -2).forward()
    NeedlemanWunch(**kwargs).forward()
    # dp = DotPlot( sys.argv[1:] )
    # p = dp.dot_plot()
    # dp.plot(p)


if __name__ == "__main__":

    main()