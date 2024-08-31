import time
import datetime

def print_header_cds(tag):
    print()
    print()
    print(f"{'IGGYPOP'.center(80)}")
    print('-' * 80)
    print(f"{'Indexed golden gate DNA assembly by PCR amplification off oligo pool templates'.center(80)}")
    print('-' * 80)
    print()
    print("    CODING SEQUENCE MODE")
    print()
    print("    **the default yaml/moclo_cds_mcu.yml configuration file:")
    print("    ")
    print("          removes common gg enzyme sites")
    print("          codon optimizes with an Arabidopsis codon table using match_codon_usage")
    print("          removes internal micro-homologies (10 bp repeats) and hairpins")
    print("          ensures GC content ≤ 0.60 over a 60 bp window size")
    print("          prevents homopolymeric runs of 6 or more of the same base")
    print("          generates MoClo compatible CDSs")
    print()
    print('    use:')
    print('       "--codon_opt use_best_codon" for best codon optimization')
    print('       "--codon_opt hybrid" for higher CAI than match_codon_usage')
    print('       "--species e_coli" for coli codon optimization')
    print('       "--mode no_mods" to skip DNA chiseling; no changes to input seq')
    print('       "--mode no_hinge" to chisel sequences')
    print("       and/or set run parameters using a YAML")
    print("")
    print("     see the yaml/ folder for other common run paramaters")
    print()
    print('-' * 80)
    print(f"{'Sean Cutler            cutler@ucr.edu'.center(80)}")
    print()
    print(f"{'University of California, Riverside'.center(80)}")
    print('-' * 80)
    print()
    print()
    time.sleep(1)


def print_header_gb(tag):
    print()
    print()
    print(f"{'IGGYPOP'.center(80)}")
    print('-' * 80)
    print(f"{'Indexed golden gate DNA assembly by PCR amplification off oligo pool templates'.center(80)}")
    print('-' * 80)
    print()
    print("    GENBANK FILE MODE")
    print()
    print("    Sequence constraints and optimizations are set as annotations in your")
    print("    genbank file; these annotations can be added automatically using the")
    print("    'iggypop.py format' option, which uses a yaml file to set run parameters.")
    print()
    print("    Parameters for hinging sequences can also be set as command line options.")
    print()
    print('-' * 80)
    print(f"{'Sean Cutler            cutler@ucr.edu'.center(80)}")
    print()
    print(f"{'University of California, Riverside'.center(80)}")
    print('-' * 80)
    print()
    print()
    time.sleep(1)
    return


def print_header_format(tag):
    print()
    print()
    print(f"{'IGGYPOP'.center(80)}")
    print('-' * 80)
    print("indexed golden gate assembly of fragmented sequences using PCR amplification")
    print(f"{'from oligo pool templates'.center(80)}")
    print('-' * 80)
    print()
    print("    GENBANK FILE FORMATTING MODE ")
    print()
    print("    Sequence constraints and optimizations are set as annotations in your")
    print("    input genbank file; these annotations can be added automatically using")
    print("    the 'format' option.")
    print()
    print('-' * 80)
    print(f"{'Sean Cutler            cutler@ucr.edu'.center(80)}")
    print()
    print(f"{'University of California, Riverside'.center(80)}")
    print('-' * 80)
    print()
    time.sleep(1)
    return


def nerd_alert():
    print()
    print()
    print("                            P7  .:!7!~.    .~?!:   .~:     ")
    print("                            J?7?7~:.    .^77~::^:    ^~    ")
    print("                        .^!?7!~.     .~!!~.     :^:   :~   ")
    print("                    .^!7!~:     ~5!!??77!~!!!!!77!J??7 ^~  ")
    print("                 :!??~:      :!?JY~:.:.:::::::::^:::~7. J. ")
    print("              :!7!^     .:~7??!:                      ~~5. ")
    print("           :!?!:    .^!77!^.Y^                        .J7. ")
    print("        .~7!^    :~77~:     ?^                       .J:?. ")
    print("      :!!^   :?!?YJ7~~~!!!!!~!~^!~~!~~~~!~~~~!~?!!!!?P^:J  ")
    print("    :!~.   :!??^ .......::::~:...::::::::::::^:^^~YJ^: J~  ")
    print("  .!~.   ^77^               Y:                 .77:   ??   ")
    print(" .?.   ^!~.                 P^              .~7!:   .?7    ")
    print(" ?:::~Y~:..::::::::.........~:............^??!.   .!J^     ")
    print(":7 7J7!~~~~~!!!!~!!~~~~~!!!!~!!!!!!!~!7?Y?7P^   :7?~       ")
    print("!^~!                        J:    .^!??~.    .~??^         ")
    print("!!7                         5: .^!7!^.     ^7?7:           ")
    print("??~                        .?!!!~:      :!?7^              ")
    print("J~~!..:::::....:.......:^77~^:      .^!7~:                 ")
    print(":?!!7?!!777!!!!7!!7?J??~~Y:      :~77~:                    ")
    print(" ~^  :~^.      .^!?7^.       :^!!!^.                       ")
    print("  ~^   .^^^. ^77~:      .:~!?!^.                           ")
    print("   :^.    ^?!~.      :~!!~^:5.                             ")
    print("    .^:^775^     .^!!~:.   .P.                             ")
    print("     .!7^ :   :~7!~~:      :G.                             ")
    print("    ~~.    .~7!^.   ::^^.  ^B.                             ")
    print("  :7:    ^7!^^^.       .^~^:G.                             ")
    print()
    print()
    return


def generate_output_filename(stem):
    """
    Generate an output filename with a timestamp.

    Parameters:
    stem (str): The base name for the output file.

    Returns:
    ofile (str): The output filename.
    tag (str): The tag for the output.
    """
    tag = f"{stem}_{datetime.datetime.now().strftime('%H%M%S')}"
    ofile = f"out/{tag}"
    return ofile, tag
