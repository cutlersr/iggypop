"""
This script handles initialization and configuration for the sequence optimization pipeline, 

`initialize`: Manages the setup process by parsing command-line arguments, loading default 
values based on the run type (CDS or GenBank mode), and creating output directories 
with unique tags to avoid overwriting existing data.

`parse_arguments`: Defines and parses command-line arguments for the sequence optimization 
pipeline, allowing for customization of input/output paths, codon optimization preferences, 
assembly options, and miscellaneous settings.

`set_defaults`: Sets default values for various parameters based on the run type, providing 
baseline settings for codon optimization, overhang sequences, assembly options, and other 
pipeline configurations.

`load_config_and_set_globals`: Loads YAML configuration, merges it with default settings, 
and applies any command-line overrides, providing a final set of parameters for the run.

`copy_analysis_assets`: Copies the specified input FASTA and YAML files, along with analysis 
scripts, to the output directory to ensure all resources are available for downstream analysis.

Additional Functions:
- `read_log_and_identify_failures`: Reads the log file to identify any sequences that 
  failed during processing.
- `check_ext_overhangs`: Verifies whether external overhangs meet fidelity requirements, 
  warning the user if low-fidelity overhangs are requested.
"""

import os
import sys
import argparse
import yaml
import shutil

def initialize(run_type="cds"):
    import sys  # Ensure sys is imported if not already

    # Parse arguments based on the run_type and set defaults
    if run_type == "cds":
        args = parse_arguments(run_type, default_yml='yaml/domesticate_cds.yml')
    elif run_type == "gb" or run_type == "format":
        args = parse_arguments(run_type, default_yml='yaml/domesticate_gb.yml')
    else:
        print("Please use a valid run type (gb or cds)")
        sys.exit(1)

    # Load default values
    default_values = set_defaults(run_type)

    # Load YAML configuration and apply command-line arguments
    updated_args, tag = load_config_and_set_globals(args, default_values)

    if updated_args.i:
        fasta_file_path = f'{updated_args.i}'
    else:
        sys.exit("Please provide a valid input file.")

    # **Add this check for input file existence**
    if not os.path.isfile(fasta_file_path):
        sys.exit(f"Error: Input file '{fasta_file_path}' does not exist.")

    # Generate output directory and log file
    tag = f"{updated_args.o}"
    base_dir = 'out'
    new_tag = tag

    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    counter = 1
    while os.path.exists(os.path.join(base_dir, new_tag)):
        new_tag = f'{tag}_{counter}'
        counter += 1
    tag = new_tag

    if updated_args.taxIDs:
        print_taxIDs()
        sys.exit(0)

    # Copy input files and analysis scripts to the results folder
    yml_file_path = f'{updated_args.yml}'
    script_paths = [
        "iggypop.py", "iggypop/iggypop_cds.py", "iggypop/chisel_hinge.py",
        "iggypop/pop_helpers.py", "iggypop/intron_stuff.py"
    ]

    new_dir = os.path.join(base_dir, tag)
    ofile = os.path.join(new_dir, tag)
    log_file_path = os.path.join(new_dir, f"log.txt")

    if run_type != "format":
        os.makedirs(new_dir)
        copy_analysis_assets(fasta_file_path, yml_file_path, tag, script_paths)

    return tag, ofile, log_file_path, updated_args



def parse_arguments(run_type="cds", default_yml='yaml/domesticate_cds.yml'):

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
        help=' (str) Path to the file containing the\n'
             ' sequences to pop\n'
    )

    general_group.add_argument(
        '--o', type=str, default="pops", metavar='',
        help=' (str) Output file stem; a timestamp is\n'
             ' added after the stem name. Output data is\n'
             ' saved to the out/ folder.\n'
    )

    general_group.add_argument(
        '--yml', default=default_yml, type=str, metavar='',
        help=' (str) YAML file that defines the run paramaters\n'
             ' Parameters set on the command line override those in\n'
             ' the yaml file. Default is yaml/domesticate_cds.yml.\n'
             ' See the yaml/ folder for many common workflows.\n'
    )

    general_group.add_argument(
        '--mode', type=str, choices=['chisel', 'no_mods', 'no_hinge'], metavar='',
        help=' (str) Mode of operation:\n'
             '  chisel - Default mode, modifies the input\n'
             '  sequence\n'
             '  no_mods - No modifications, just hinges\n'
             '  and barcodes\n'
             '  no_hinge - Chisels the sequences without\n'
             '  hinging\n'
    )

    general_group.add_argument(
        '--reports', action='store_true', default=True,
        help=' (bool) Enable DNAchisel reports\n'
    )

    general_group.add_argument(
        '--no-reports', action='store_false', dest='reports',
        help='Disable DNAchisel reports\n'
    )


    # Codon optimization arguments
    if run_type == 'cds' or run_type == 'format':
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
            help=' (float) Percent sequence divergence from\n'
                 ' the perfect use_best_codon sequence to\n'
                 ' target in hybrid codon optimization mode.\n'
                 ' Default is 20 pct.\n'
        )

        codon_group.add_argument(
            '--species', type=str, metavar='',
            help=' (str) Specify the codon table to use for\n'
                 ' optimization; default is arabidopsis.\n'
                 ' Can use NCBI taxIDs for other species.\n'
                 ' Connects to Kazusa for everything except\n'
                 ' e_coli, s_cerevisiae, b_subtilis, and a\n'
                 ' few other model systems.\n'
        )

        codon_group.add_argument(
            '--codon_tbl', type=str, choices=[
                'cocoputs', 'kazusa'
            ], metavar='',
            help='(str) Codon table lookup:\n'
                 ' Use cocoputs or kazusa codon data. You can\n'
                 ' use short names like e_coli, a_thaliana,\n'
                 ' b_subtilis, etc. in the --species parameter;\n'
                 ' you can also provide a taxID; cocoputs does not\n'
                 ' require an internet connection and uses\n'
                 ' more recent and extensive codon data than\n'
                 ' Kazusa.\n'
        )

        codon_group.add_argument(
            '--original_species', type=str, default="none", metavar='',
            help=' (str) Source codon table to harmonize\n'
                 ' with, used for harmonize_rca.\n'
        )

        codon_group.add_argument(
            '--deintronize', type=str, choices=[
                'on', 'off'
            ], metavar='',
            help='(str) Codon table lookup:\n'
                 ' enable deinitronization mode; experimental fature\n'
        )

        codon_group.add_argument(
            '--require_orf', type=str, choices=[
                'on', 'off'
            ], metavar='',
            help='(str) ignore the requirements of an orf\n'
        )

    # Codon optimization arguments
    tweak_group = parser.add_argument_group('Repeated Design Options')

    tweak_group.add_argument(
        '--repeats', type=int, default=1, metavar='',
        help=' (int) Number of chiseled sequences to\n'
             ' create per input sequence. Default is 1.\n'
             ' Useful if you want to optimize by tweaking.\n'
    )

    # Assembly options
    assembly_group = parser.add_argument_group('Assembly Options')


    assembly_group.add_argument(
        '--two_step_length', type=int, metavar='',
            help='(int) maximum length for two-step fragments.\n'
            ' remember to use a two-step yaml\n'
    )

    assembly_group.add_argument(
        '--oligo_length', type=int, metavar='',
        help=' (int) Maximum length of oligo allowed\n'
             ' Default is 250'
    )

    assembly_group.add_argument(
        '--primer_length', type=int, metavar='',
        help=' (int) Barcode primer length\n'
             ' Default is 18'
    )


    assembly_group.add_argument(
        '--ext_overhangs', type=str, nargs='+', metavar='',
        help=' (list) External overhangs used for cloning.\n'
             ' Defaults are AATG and GCTT. This parameter\n'
             ' is used to prevent the external cloning\n'
             ' overhangs from occurring inside the internal\n'
             ' assembly junctions; you need to separately\n'
             ' define corresponding end sequences using\n'
             ' base_5p_end and base_3p_end.\n'
    )

    assembly_group.add_argument(
        '--base_5p_end', type=str, metavar='',
        help=' (str) Sequence to append to the 5\' end of\n'
             ' the chiseled CDS; default is\n'
             ' "AATGCGGTCTCTA" for MoClo compatible CDSs.\n'
    )

    assembly_group.add_argument(
        '--base_3p_end', type=str, metavar='',
        help=' (str) Sequence to append to the 3\' end of\n'
             ' the chiseled CDS; default is\n'
             ' "GCTTAGAGACCGCTT" for MoClo compatible CDSs.\n'
    )

    assembly_group.add_argument(
        '--pcr_5p_cut', type=str, metavar='',
        help=' (str) Sequence appended to the beginning\n'
             ' of each oligo for gg cloning. Default is\n'
             ' BsmBI --> "CGTCTCA".\n'
    )

    assembly_group.add_argument(
        '--pcr_3p_cut', type=str, metavar='',
        help=' (str) Sequence appended to the end of each\n'
             ' oligo for gg cloning. Default is BsmBI -->\n'
             ' "AGAGACG".\n'
    )

    assembly_group.add_argument(
        '--max_fragments', type=int, metavar='',
        help=' (int) Max number of fragments per PCR\n'
             ' - default 18\n'
    )

    assembly_group.add_argument(
        '--primer_index', type=int, metavar='',
        help=' (int) Where to start adding primers from\n'
             ' the indexing primer set\n'
    )

    assembly_group.add_argument(
        '--n_tries', type=int, metavar='',
        help=' (int) Number of potential overhang sets to\n'
             ' consider for fragment junctions. Default is\n'
             ' 50. The highest fidelity solution obtained\n'
             ' from n_tries solutions is used for assemblies.\n'
             ' Decreasing will speed things up with a\n'
             ' relatively minor impact on assembly fidelity\n'
             ' scores.\n'
    )

    assembly_group.add_argument(
        '--radius', type=int, metavar='',
        help=' (int) Allowable distance from ideal cut\n'
             ' sites for selecting OHs. Default 8; \n' 
             ' decreasing this value will increase\n'
             ' the average oligo size but may lower\n'
             ' the fidelity of assemblies slightly.\n'
             ' Increasing will reduce average oligo length \n'
             ' but provide higher ligation fidelities.\n'
    )

    # Miscellaneous options
    misc_group = parser.add_argument_group('Miscellaneous Options')

    misc_group.add_argument(
        '--seed', dest='seed', type=int, default=False, metavar='',
        help=' (int) Seed value for random number generator.\n'
             ' Set this to reproduce an older run\'s results.\n'
             ' Default is False.\n'
    )

    misc_group.add_argument(
        '--index_primers', type=str, metavar='',
        help=' (str) Path to the index primers file. Default\n'
             ' is "data/10K_primers_renamed.csv". Change to use\n'
             ' custom barcodes'
    )

    misc_group.add_argument(
        '--fidelity_data', type=str, metavar='',
        help=' (str) Path to OH fidelity data from Potapov et al.; \n'
             ' default = 18 hrs 25 oC T4 ligase data;\n'
             ' see the data folder for other options.\n'
    )

    misc_group.add_argument(
        '--ohsets', type=str, metavar='',
        help=' (str) High-fidelity overhang sets used for\n'
             ' generating high-fidelity assemblies.\n'
    )

    # Quiet/verbose flags: quiet defaults to “on”, --verbose sets it to “off”
    misc_group.add_argument(
        '--quiet',
        choices=['on', 'off'],
        default='on',
        metavar='',
        help='Minimize prints to terminal (default: on)'
    )
    misc_group.add_argument(
        '--verbose',
        action='store_const',
        const='off',
        dest='quiet',
        help='Turn off quiet mode (i.e. verbose output)'
    )
    
    misc_group.add_argument(
        '--taxIDs', action='store_true',
        help=' (bool) Print the list of common organisms and\n'
             ' their NCBI taxIDs.\n'
    )

    return parser.parse_args()


def set_defaults(run_type):
    default_values = {
        'base_5p_end': 'AATGCGGTCTCTA',
        'base_3p_end': 'GCTTAGAGACCGCTT',
        'ext_overhangs': ['AATG', 'GCTT'],
        'allowed_chars': "ATGC",
        'two_step': "off",
        'two_step_length': 1200,
        'two_step_5p_end': 'AATGCGTCTCA',
        'two_step_3p_end': 'AGAGACGGCTT',
        'pcr_5p_cut': "CGTCTCA",    
        'pcr_3p_cut': "AGAGACG",  
        'primer_length': 18,
        'oligo_length': 250,
        'require_orf': 'on',         
        'fidelity_data': 'data/FileS03_T4_18h_25C.xlsx',
        'ohsets': 'data/hingesets.xlsx',
        'mode': 'chisel',
        'species': 'arabidopsis',
        'primer_index': 1,
        'codon_tbl': "cocoputs",
        'original_species': 'none',
        'n_tries': 10,
        'max_fragments': 50,
        'radius': 8,
        'index_primers': 'data/indexsets.csv',
        'deintronize': 'off',
        'taxIDs': False,
        'seed': False,
        'pct': 20,
        'repeats': 1,
        'codon_opt': 'none',
        'quiet': 'on',
        'reports': "True"
    }

    if run_type == "gb":
        default_values['species'] = 'e_coli'
        default_values['codon_tbl'] = 'kazusa'

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
    if 'yml' in args_dict and args_dict['yml']:
        yml_path = args_dict['yml']
        if os.path.exists(yml_path):
            with open(yml_path, 'r') as ymlfile:
                try:
                    config = yaml.safe_load(ymlfile) or {}
                except yaml.YAMLError as e:
                    sys.exit(f"Error: Failed to parse YAML file '{yml_path}'.\n{e}")
        else:
            sys.exit(f"Error: YAML configuration file '{yml_path}' does not exist.")
    else:
        sys.exit("Error: No YAML configuration file specified.")

    # Merge values: start with defaults, then YAML config, then command-line args
    updated_values = {**default_values, **config}

    # Overwrite only with non-None values from command-line args
    for key, value in args_dict.items():
        # For store_true/store_false type args, only override if explicitly provided
        if isinstance(value, bool) and value is None:
            continue
        if value is not None:
            updated_values[key] = value

    # Check for missing required parameters
    missing_params = [key for key, value in updated_values.items() if value is None]
    if missing_params:
        sys.exit(f"Error: Missing required parameters: {', '.join(missing_params)}")

    # Convert updated_values to Namespace for easy access
    updated_args = argparse.Namespace(**updated_values)
    tag = updated_values.get('o', 'default_tag')

    return updated_args, tag


def copy_analysis_assets(fasta_input, yml_input, tag, script_paths):
    """
    Copy analysis assets to the target directory.

    Parameters:
    fasta_input (str): Path to the input FASTA file.
    yml_input (str): Path to the input YAML file.
    tag (str): Tag for the output directory.
    script_paths (list): List of paths to analysis scripts.
    """
    target_directory = f"out/{tag}/assets"
    os.makedirs(target_directory, exist_ok=True)

    # Copy the fasta and yml files
    shutil.copy(fasta_input, os.path.join(target_directory, os.path.basename(fasta_input)))

    if yml_input is not None:
        shutil.copy(yml_input, os.path.join(target_directory, os.path.basename(yml_input)))

    # Copy the specified analysis scripts
    for script_path in script_paths:
        shutil.copy(script_path, os.path.join(target_directory, os.path.basename(script_path)))

def load_data(fidelity_data, ohsets, ext_overhangs, fasta_file, log_file):
    # Load Potapov data for calculating overhang fidelities
    try:
        potapov_data = pd.read_excel(fidelity_data)
    except Exception as e:
        log_and_exit(
            f"FAILED to load fidelity data.\nException occurred: {e}\n",
            log_file, exit_code=2
        )

    # Load high-fidelity overhang sets
    try:
        overhang_sets = get_overhang_sets(ohsets, ext_overhangs)
    except Exception as e:
        log_and_exit(
            f"FAILED to get overhang sets.\nException occurred: {e}\n",
            log_file, exit_code=3
        )

    # Read input sequences
    sequences = read_fasta(fasta_file)
    if not sequences:
        log_and_exit(
            "No sequences found in the input file.", log_file, exit_code=4
        )

    return potapov_data, sequences, overhang_sets


def read_log_and_identify_failures(tag):
    """
    Read the log file for a given tag, identify any failed sequences, and print
    the total number of failed sequences along with their identifiers.

    Parameters:
    tag (str): The tag used to identify the specific log file in the results 
               directory.
    """
    log_file_path = f"out/{tag}/log.txt"
    failed_sequences = []

    try:
        # Attempting to open and read the log file
        with open(log_file_path, "r") as log_file:
            for line in log_file:
                if "No solutions sets found for" in line:
                    # Extracting the sequence ID from the log message
                    seq_id = line.split("No solutions sets found for ")[1].strip()
                    failed_sequences.append(seq_id)

        # Printing the total number of failed sequences and their identifiers
        print('*' * 80,)
        if failed_sequences:
            print("\nFailed sequences:", ", ".join(failed_sequences), "\n")
        else:
            print("\nNo failed sequences")

    except FileNotFoundError:
        print(f"Log file not found: {log_file_path}")

def check_ext_overhangs(overhangs):
    overhang_set = {'GGAC', 'CGCG', 'GGCC', 'GCCA', 'CCGC', 'GGCG', 'GTCG', 'TCGC', 'GGGT', 'GGGG', 'TGCG', 'TAAA', 
                    'GCAT', 'GGCT', 'TTTA', 'GCGG', 'TATA', 'GGAT', 'CCCC', 'GGGC', 'GCGT', 'TTAA', 'GGAG', 'GGTG'}
    
    for overhang in overhangs:
        if overhang in overhang_set:
            print(f'')             
            print(f'The external overhangs you would like to use – {overhangs} – have low ligation fidelity.') 
            print(f'')             
            print(f"The following cloning overhangs wont't work with the default data/hingesets.xlsx sets:")
            print(f"'GGAC', 'CGCG', 'GGCC', 'GCCA', 'CCGC', 'GGCG', 'GTCG', 'TCGC', 'GGGT', 'GGGG', 'TGCG'")
            print(f"'TAAA', 'GCAT', 'GGCT', 'TTTA', 'GCGG', 'TATA', 'GGAT', 'CCCC', 'GGGC', 'GCGT', 'TTAA'")
            print(f"'GGAG', 'GGTG'")
            print(f'')             
            print(f'If you have to use {overhangs}, you can use iggypop.py gagga to design a new overhang sets,')
            print(f'or provide your own custom overhang sets in place of the default data/hf_oh_sets.xlsx.')             
            print(f'')             
            return True
    return False

