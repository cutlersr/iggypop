# initialization.py

import os
import sys
import yaml
import argparse
import pandas as pd
from pop_helpers import log_and_exit

def initialize_configuration_and_output(default_yml='yaml/moclo_cds.yml'):
    """
    Parse arguments, set defaults, load configuration, and generate output directory.

    Parameters:
    default_yml (str): Path to the default YAML file.

    Returns:
    args (Namespace): Parsed arguments.
    updated_defaults (Namespace): Updated arguments with combined values.
    tag (str): Tag for the configuration.
    ofile (str): Output file path.
    potapov_data (DataFrame): Loaded Potapov data.
    overhang_sets (dict): Loaded high-fidelity overhang sets.
    log_file (file object): Opened log file for writing.
    """
    def parse_arguments(default_yml):
        parser = argparse.ArgumentParser(
            description=(
                'Optimize, fragment, select high-fidelity overhangs, '
                'and barcode sequences for gg assemblies'
            ),
            formatter_class=argparse.RawTextHelpFormatter
        )

        # General arguments
        general_group = parser.add_argument_group('General Options')

        general_group.add_argument(
            '--i', type=str, metavar='',
            help='(str) Path to the file containing the\n'
                 'sequences to pop\n'
        )

        general_group.add_argument(
            '--o', type=str, default="pops", metavar='',
            help='(str) Output file stem; a timestamp is\n'
                 'added after the stem name. Output data is\n'
                 'saved to the out/ folder.\n'
        )

        general_group.add_argument(
            '--yml', default=default_yml, type=str, metavar='',
            help='(str) YAML file that defines the\n'
                 'restriction enzymes to avoid, codon\n'
                 'optimization, etc.; default is\n'
                 'yaml/default.yml. See the yaml/ folder\n'
                 'for many common workflows.\n'
        )

        general_group.add_argument(
            '--mode', type=str, choices=['chisel', 'no_mods', 'no_hinge'], metavar='',
            help='(str) Mode of operation:\n'
                 '  chisel - Default mode, modifies the input\n'
                 '  sequence\n'
                 '  no_mods - No modifications, just hinges\n'
                 '  and barcodes\n'
                 '  no_hinge - Chisels the sequences without\n'
                 '  hinging\n'
        )

        general_group.add_argument(
            '--reports', action='store_true', 
            help='(bool) Enable DNAchisel reports\n'
        )

        # Codon optimization arguments
        codon_group = parser.add_argument_group('Codon Optimization Options')

        codon_group.add_argument(
            '--codon_opt', type=str, choices=[
                'use_best_codon', 'match_codon_usage', 'harmonize_rca', 'hybrid', 
                'none'
            ], metavar='',
            help='(str) Codon optimization method:\n'
                 '  use_best_codon - Use the best codon for\n'
                 '  optimization\n'
                 '  match_codon_usage - Match codon usage\n'
                 '  harmonize_rca - Harmonize with RCA\n'
                 '  hybrid - Use best codon with sequence\n'
                 '  diversification\n'
                 '  none - No codon optimization\n'
        )

        codon_group.add_argument(
            '--pct', type=float, default=20, metavar='',
            help='(float) Percent sequence divergence from\n'
                 'the perfect use_best_codon sequence to\n'
                 'target in hybrid codon optimization mode.\n'
                 'Default is 20 pct.\n'
        )

        codon_group.add_argument(
            '--species', type=str, metavar='',
            help='(str) Specify the codon table to use for\n'
                 'optimization; default is arabidopsis.\n'
                 'Can use NCBI taxIDs for other species.\n'
                 'Connects to Kazusa for everything except\n'
                 'e_coli, s_cerevisiae, b_subtilis, and a\n'
                 'few other model systems.\n'
        )

        codon_group.add_argument(
            '--cocoputs', action='store_true',
            help='(bool) Use cocoputs codon data. You can\n'
                 'use short names like e_coli, a_thaliana,\n'
                 'b_subtilis, etc. in the --species parameter;\n'
                 'you can also provide a taxID; does not\n'
                 'require an internet connection and uses\n'
                 'more recent and extensive codon data than\n'
                 'Kazusa.\n'
        )

        codon_group.add_argument(
            '--original_species', type=str, default="none", metavar='',
            help='(str) Source codon table to harmonize\n'
                 'with, used for harmonize_rca.\n'
        )

        codon_group.add_argument(
            '--deintronize', action='store_true',
            help='(bool) Enable deintronization\n'
                 '(experimental feature)\n'
        )

        codon_group.add_argument(
            '--repeats', type=int, default=1, metavar='',
            help='(int) Number of chiseled sequences to\n'
                 'create per input sequence. Default is 1.\n'
                 'Useful if you want to optimize by tweaking.\n'
        )

        codon_group.add_argument(
            '--tweak_cai', type=float, default=0.8, metavar='',
            help='(float) this sets a minimum CAI value to\n'
                 'use when you are optimiing by tweaking.\n'
                 'default = 0.8. *** no_hinge mode only.\n'
        )

        codon_group.add_argument(
            '--tweak_n', type=int, default=False, metavar='',
            help='(int) Number of tweaked sequences you\n'
                 'desire. A typical protocol is to run.\n'
                 '25 repeats and keep the 5 most divergent\n' 
                 'sequences with CAI values greater than 0.8\n'
                 ' *** no_hinge mode only\n'
        )

        # Assembly options
        assembly_group = parser.add_argument_group('Assembly Options')

        assembly_group.add_argument(
            '--two_step', action='store_true', 
            help='(bool) Enable two-step assemblies; if you\n'
                 'use this mode, use the two_step.yml file\n'
                 'and pay close attention to ends & enzymes.\n'
        )

        assembly_group.add_argument(
            '--max_fragments', type=int, metavar='',
            help='(int) Max number of fragments per PCR\n'
                 '- default 6\n'
        )

        assembly_group.add_argument(
            '--segment_length', type=int, metavar='',
            help='(int) Maximum length of segments to be\n'
                 'created from input sequence; 200 for\n'
                 'default to allow for barcodes & cut sites,\n'
                 'assuming 250 bp oligos.\n'
        )

        assembly_group.add_argument(
            '--ext_overhangs', type=str, nargs='+', metavar='',
            help='(list) External overhangs used for cloning.\n'
                 'Defaults are AATG and GCTT. This parameter\n'
                 'is used to prevent the external cloning\n'
                 'overhangs from occurring inside the internal\n'
                 'assembly junctions; you need to separately\n'
                 'define corresponding end sequences using\n'
                 'base_5p_end and base_3p_end.\n'
        )

        assembly_group.add_argument(
            '--base_5p_end', type=str, metavar='',
            help='(str) Sequence to append to the 5\' end of\n'
                 'the chiseled CDS; default is\n'
                 '"AATGCGGTCTCTA" for MoClo compatible CDSs.\n'
        )

        assembly_group.add_argument(
            '--base_3p_end', type=str, metavar='',
            help='(str) Sequence to append to the 3\' end of\n'
                 'the chiseled CDS; default is\n'
                 '"GCTTAGAGACCGCTT" for MoClo compatible CDSs.\n'
        )

        assembly_group.add_argument(
            '--pcr_5p_cut', type=str, metavar='',
            help='(str) Sequence appended to the beginning\n'
                 'of each oligo for gg cloning. Default is\n'
                 'BsmBI --> "CGTCTCA".\n'
        )

        assembly_group.add_argument(
            '--pcr_3p_cut', type=str, metavar='',
            help='(str) Sequence appended to the end of each\n'
                 'oligo for gg cloning. Default is BsmBI -->\n'
                 '"AGAGACG".\n'
        )

        assembly_group.add_argument(
            '--primer_index', type=int, metavar='',
            help='(int) Where to start adding primers from\n'
                 'the indexing primer set\n'
        )

        assembly_group.add_argument(
            '--n_tries', type=int, metavar='',
            help='(int) Number of potential overhang sets to\n'
                 'consider for fragment junctions. Default is\n'
                 '50. The highest fidelity solution obtained\n'
                 'from n_tries solutions is used for assemblies.\n'
                 'Decreasing will speed things up with a\n'
                 'relatively minor impact on assembly fidelity\n'
                 'scores.\n'
        )

        assembly_group.add_argument(
            '--radius', type=int, metavar='',
            help='(int) Allowable distance from ideal cut\n'
                 'sites. Default 8; decreasing this value will\n'
                 'increase the average oligo size but may lower\n'
                 'the fidelity of assemblies slightly.\n'
        )

        # Miscellaneous options
        misc_group = parser.add_argument_group('Miscellaneous Options')

        misc_group.add_argument(
            '--seed', dest='seed', type=int, default=False, metavar='',
            help='(int) Seed value for random number generator.\n'
                 'Set this to reproduce an older run\'s results.\n'
                 'Default is False.\n'
        )

        misc_group.add_argument(
            '--index_primers', type=str, metavar='',
            help='(str) Path to the index primers file. Default\n'
                 'is "data/10K_primers_renamed.csv". Change to use\n'
                 'custom barcodes; our defaults are built around\n'
                 '18 bp primers; update segment length as needed\n'
                 'if the barcode primer length changes.\n'
        )

        misc_group.add_argument(
            '--fidelity_data', type=str, metavar='',
            help='(str) OH fidelity data from Potapov et al.; see\n'
                 'the data folder for other options.\n'
        )

        misc_group.add_argument(
            '--ohsets', type=str, metavar='',
            help='(str) High-fidelity overhang sets used for\n'
                 'generating high-fidelity assemblies.\n'
        )

        misc_group.add_argument(
            '--taxIDs', action='store_true',
            help='(bool) Print the list of common organisms and\n'
                 'their NCBI taxIDs.\n'
        )

        return parser.parse_args()

    def set_defaults():
        default_values = {
            'base_5p_end': 'AATGCGGTCTCTA',
            'base_3p_end': 'GCTTAGAGACCGCTT',
            'ext_overhangs': ['AATG', 'GCTT'],
            'allowed_chars': "ATGC",
            'two_step': False,
            'two_step_length': 1104,
            'two_step_5p_end': 'AATGGGTCTCA',
            'two_step_3p_end': 'TGAGACCGCTT',
            'segment_length': 200,
            'fidelity_data': 'data/FileS03_T4_18h_25C.xlsx',
            'ohsets': 'data/hf_oh_sets.xlsx',
            'mode': 'chisel',
            'species': 'arabidopsis',
            'primer_index': 1,
            'cocoputs': False,
            'original_species': False,
            'n_tries': 150,
            'max_fragments': 7,
            'radius': 8,
            'index_primers': 'data/10K_primers_renamed.csv',
            'pcr_5p_cut': "CGTCTCA",
            'pcr_3p_cut': "AGAGACG",
            'cutoff': 0.45,
            'tweak_n': False,
            'tweak_cai': 0.8,
            'deintronize': False,
            'original_species': 'none',
            'taxIDs': False,
            'seed': False,
            'pct': 20,
            'repeats': 1,
            'codon_opt': 'match_codon_usage'
        }
        return default_values

    def load_config_and_set_globals(args, default_values):
        """
        Load the YAML configuration and update global settings.

        Parameters:
        args (Namespace or dict): Arguments from argparse or a dictionary.
        default_values (dict): Default values for the configuration.

        Returns:
        updated_args (Namespace): Updated arguments with combined values.
        tag (str): Tag for the configuration.
        """
        # Convert args to a dictionary if it's a Namespace from argparse
        args_dict = vars(args) if isinstance(args, argparse.Namespace) else args

        # Load the YAML configuration
        config = {}
        if 'yml' in args_dict and args_dict['yml'] and os.path.exists(args_dict['yml']):
            with open(args_dict['yml'], 'r') as ymlfile:
                config = yaml.safe_load(ymlfile) or {}

        updated_values = {**default_values, **config}

        for key, value in args_dict.items():
            # Update the value only if it's not None, to avoid overwriting with None
            if value is not None:
                updated_values[key] = value

        missing_params = [key for key, value in updated_values.items() if value is None]
        if missing_params:
            raise ValueError(f"Missing required parameters: {', '.join(missing_params)}")

        updated_args = argparse.Namespace(**updated_values)
        tag = updated_values.get('tag', 'default_tag')

        return updated_args, tag



    args = parse_arguments(default_yml)
    default_values = set_defaults()
    updated_defaults, tag = load_config_and_set_globals(args, default_values)

    # Update globals with the updated default values
    globals().update(vars(updated_defaults))

    # Generate output directory and log file
    tag = f"{args.o}"
    base_dir = 'out'
    new_tag = tag

    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    counter = 1
    while os.path.exists(os.path.join(base_dir, new_tag)):
        new_tag = f'{tag}_{counter}'
        counter += 1
    tag = new_tag

    new_dir = os.path.join(base_dir, tag)
    os.makedirs(new_dir)
    ofile = os.path.join(new_dir, tag)

    # Open log file
    log_file = open(f"{new_dir}/{tag}_log.txt", "a")

    # Load Potapov data for calculating overhang fidelities
    try:
        potapov_data = pd.read_excel(updated_defaults.fidelity_data)
    except Exception as e:
        log_and_exit(
            f"FAILED to load fidelity data.\nException occurred: {e}\n",
            log_file, exit_code=2
        )

    # Load high-fidelity overhang sets
    try:
        overhang_sets = get_overhang_sets(updated_defaults.ohsets, updated_defaults.ext_overhangs)
    except Exception as e:
        log_and_exit(
            f"FAILED to get overhang sets.\nException occurred: {e}\n",
            log_file, exit_code=3
        )

    return args, updated_defaults, tag, ofile, potapov_data, overhang_sets, log_file
