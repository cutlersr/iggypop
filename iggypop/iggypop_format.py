import argparse
import yaml
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import warnings
from initialization import *
from pop_helpers import *

def annotate_genbank():
    try:
        # Load the GenBank file
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(i, "genbank"))  # 'i' needs definition where it's used or passed
            for warning in w:
                custom_warning_handler(warning.message, warning.category, warning.filename, warning.lineno)
        if not records:
            raise ValueError("No records found in the GenBank file.")

        # Load the YAML file
        with open(yml, 'r') as yf:  # 'yaml' here is ambiguous as it is also the name of the imported module
            yaml_data = yaml.safe_load(yf)

        # Annotate records using the parsed YAML data
        annotated_records = annotate_from_yaml(yaml_data, species, codon_opt, records, pct)  # 'pct' also needs definition

        # Write the modified records back to the output file
        with open(o, "w") as output_handle:
            SeqIO.write(annotated_records, output_handle, "genbank")

        # Post-process the output file to replace the specific label
        with open(o, "r") as file:
            content = file.read()

        content = content.replace(
            f'                     /label="~EnforceChanges(minimum_percent={pct})"',
            f'                     /label="@EnforceChanges(minimum={pct}%)"'
        )

        # Write the modified content back to the file
        with open(o, "w") as file:
            file.write(content)

        print(f"Annotation completed successfully. Output written to {o}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print("Please check the inputs and try again.")

def annotate_from_yaml(yaml_data, species, codon_opt, records, pct):
    regulatory_features = [
        "promoter", "terminator", "RBS", "regulatory",
        "rep_origin", "protein_bind", "tRNA", "oriT",
        "stem_loop", "enhancer", "5'UTR", "3'UTR",
        "polyA_signal", "LTR", "misc_binding", "tRNA",
        "rRNA", "misc_RNA"
    ]

    for record in records:
        new_features = []  # Collect new features here

        # Handle regulatory features
        for feature in record.features:
            if feature.type in regulatory_features:
                avoid_changes_feature = SeqFeature(
                    feature.location,
                    type="misc_feature",
                    qualifiers={"label": "@AvoidChanges"}
                )
                new_features.append(avoid_changes_feature)

            elif feature.type in ["CDS", "ORF"]:
                misc_feature = SeqFeature(
                    feature.location,
                    type="misc_feature",
                    qualifiers={"label": "@EnforceTranslation"}
                )
                new_features.append(misc_feature)

                if codon_opt == 'match_codon_usage':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~match_codon_usage(species={species})"}
                    )
                    new_features.append(codon_feature)

                elif codon_opt == 'use_best_codon':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~use_best_codon(species={species})"}
                    )
                    new_features.append(codon_feature)

                elif codon_opt == 'hybrid':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~use_best_codon(species={species})"}
                    )
                    new_features.append(codon_feature)
                    enforce_changes_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~EnforceChanges(minimum_percent={pct})"}
                    )
                    new_features.append(enforce_changes_feature)

                elif codon_opt == 'harmonize_rca':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~match_codon_usage(species={species})"}
                    )
                    new_features.append(codon_feature)

        # Create features for constraints from YAML
        for constraint in yaml_data.get("constraints", []):
            if constraint['type'] == "EnforceTranslation":  # Skip if the constraint is of type "EnforceTranslation"
                continue
            label = create_feature_label("@", constraint)
            constraint_feature = SeqFeature(
                FeatureLocation(0, len(record)),
                type="misc_feature",
                qualifiers={"label": label}
            )
            new_features.append(constraint_feature)

        # Create features for optimizations from YAML
        for optimization in yaml_data.get("optimizations", []):
            label = create_feature_label("~", optimization)
            optimization_feature = SeqFeature(
                FeatureLocation(0, len(record)),
                type="misc_feature",
                qualifiers={"label": label}
            )
            new_features.append(optimization_feature)

        record.features.extend(new_features)
    return records


if __name__ == "__main__":
    tag, ofile, log_file_path, updated_defaults = initialize("format")
    globals().update(vars(updated_defaults))

    # Convert species to NCBI taxID for dnachisel
    if codon_tbl == "kazusa": 
        if species == "arabidopsis" or species=="a_thaliana":
            species = '3702'
            print("using TaxID 3702 for Arabidopsis thaliana")

    if codon_tbl == "cocoputs": 
        print("you can't use cocoputs with gb mode; run with `--codon_tbl kazusa` and supply a TaxID.")
        sys.exit()


    annotate_genbank()
