
"""
Uses Spliceator models to identify potential intron donor and acceptor sites in a sequence.
Adapted from the Sppliceator project: https://git.unistra.fr/nscalzitti/spliceator; original 
author  Nicolas Scalzitti
"""

import tensorflow as tf
import os
import sys
import argparse
from os.path import basename
from pathlib import Path
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from tensorflow.keras.models import load_model
from pkg_resources import resource_filename
import logging
import absl.logging
from contextlib import redirect_stdout, redirect_stderr

# Process datasets
def intron_seqs_to_avoid(data_set):
    intron_seqs = []
    intron_locs_left = []
    intron_locs_right = []
    for location, sequence_data in data_set.items():
        sequence = sequence_data[1]

        value = sequence[9:13]
        intron_seqs.append(value)

        location_int = int(location)
        intron_locs_left.append(location_int-3)
        intron_locs_right.append(location_int+3)

    return intron_seqs, intron_locs_left, intron_locs_right

def check_intron(sequence):

    logging.getLogger('tensorflow').disabled = False
    absl.logging.set_verbosity(absl.logging.ERROR)
    model_path = "data/"
    model_donor_size = 200
    model_acceptor_size = 200
    reliability_acc = 98
    reliability_don = 98
    intron = False
    all_results = []

    model_donor = tf.keras.models.load_model(model_path + f"donor_{model_donor_size}.h5")
    model_acceptor = tf.keras.models.load_model(model_path + f"acceptor_{model_acceptor_size}.h5")
    with open(os.devnull, 'w') as f:
        with redirect_stdout(f), redirect_stderr(f):
            donor_dict, acceptor_dict = evaluate(sequence, model_donor, model_acceptor, all_results, threshold=reliability_don, window_size=model_donor_size)
    
    donor_count = 0
    acceptor_count = 0

    for k, v in donor_dict.items():
        donor_count += 1

    for k, v in acceptor_dict.items():
        acceptor_count += 1

    if (donor_count + acceptor_count) > 0: 
        intron = True

    return intron, donor_dict, acceptor_dict

def one_hot_encoding(sequence):
    
    sequence = sequence.upper()
    encoded_sequence = ""
    for nuc in sequence:
        if nuc == "A":
            encoded_sequence += "1    0    0    0    "
        elif nuc == "C":
            encoded_sequence += "0    1    0    0    "
        elif nuc == "G":
            encoded_sequence += "0    0    1    0    "
        elif nuc == "T":
            encoded_sequence += "0    0    0    1    "
        elif nuc == "N":
            encoded_sequence += "0    0    0    0    "
        else:
            encoded_sequence = ""
            break
    return encoded_sequence

def find_seq(sequence, pos, size=400):
    sequence = sequence.upper()
    seq_length = len(sequence)
    for i in range(seq_length):
        if i == pos:
            if size == 20:
                window = sequence[i:i+20]
            elif size == 80:
                window = sequence[i+30:i+50]
            elif size == 140:
                window = sequence[i+60:i+80]
            elif size == 200:
                window = sequence[i+90:i+110]
            elif size == 400:
                window = sequence[i+190:i+210]
            return window

# Returns the positions of the donor and acceptor sites
def evaluate(sequence, model_donor, model_acceptor, all_results, threshold=95, window_size=400):
    # all predictions
    donor_dict = {}
    acceptor_dict = {}
    if all_results:
        threshold = 0
    else:
        threshold = float(threshold/100)
    # sequence to analyze
    sequence = "N"*int(window_size/2) + sequence.upper() + "N"*int(window_size/2)
    # all windowed_sequences from the raw sequences
    output_donor = []
    output_acceptor = []
    d_window_size = int(window_size/4)
    # analyze of all sequence
    for i in tqdm(range(len(sequence))):
        # windowed_sequence's size is the same size as seq_train/seq_test
        windowed_sequence = sequence[i: window_size+i]
        # Encoding the target windowed sequence in one-hot
        encoded_windowed_sequence = one_hot_encoding(windowed_sequence).replace(" ","")
        # Check the length of sequences, only sequences with length=window_size are analyzed
        if len(encoded_windowed_sequence) != window_size*4:
            pass
        else:
            # Conversion
            to_add = np.array(list(encoded_windowed_sequence), dtype=int)
            # Reshaping 
            to_add2 = to_add.reshape(-1, int(window_size), 4)
            # Run the prediction on all windowed_sequences
            output_donor.append(model_donor.predict(to_add2))
            output_acceptor.append(model_acceptor.predict(to_add2))
    # remove last N
    del output_donor[-1]
    del output_acceptor[-1]
    # Donor
    for position, proba in enumerate(output_donor):
        proba = [proba[0][0], proba[0][1]]
        if all_results:
            s = find_seq(sequence, position, window_size)
            donor_dict[str(position)] = [proba, s]
        else:
            # get predicted SS
            if proba[1] > proba[0]:
                if proba[1] >= float(threshold):
                    s = find_seq(sequence, position, window_size)
                    donor_dict[str(position)] = [proba, s]
    # Acceptor
    for position, proba in enumerate(output_acceptor):
        proba = [proba[0][0], proba[0][1]]
        if all_results:
            s = find_seq(sequence, position, window_size)
            acceptor_dict[str(position)] = [proba, s]
        else:
            # get predicted SS
            if proba[1] > proba[0]:
                if proba[1] >= float(threshold):
                    s = find_seq(sequence, position, window_size)
                    acceptor_dict[str(position)] = [proba, s]
    return donor_dict, acceptor_dict
