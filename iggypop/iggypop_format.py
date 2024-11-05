"""
This script annotates GenBank files based on a default or user-provided YAML configuration, 
adding protection labels for regulatory features and applying custom constraints or optimizations.

`create_feature_label` generates feature labels for annotations based on required dnachisel prefixes 
('@' for constraints or '~' for optimizations).

`annotate_from_yaml` is responsible for reading the YAML and applying constraints and
optimizations. @AvoidChanges tags are added to regulatory features; coding sequences and ORFS
are annotated with @EnforceTranslation and any codon optimization methods if requested.

`annotate_genbank` manages the GenBank and YAML file loading, applies annotations,and writes updated 
records.

"""

import argparse
import yaml
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import warnings
from initialization import *
from pop_helpers import *

def create_feature_label(prefix, data):
    try:
        label = f"{prefix}{data['type']}"
        additional_info = [
            f"{key}={value}" for key, value in data.items() if key != 'type'
        ]
        if additional_info:
            label += f"({', '.join(additional_info)})"
        return label
    except KeyError as e:
        print(f"Missing key in data: {e}")
        print("Please ensure the YAML data contains the correct keys.")
        return f"{prefix}unknown"

def annotate_genbank(comment):
    try:
        # Load the GenBank file
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(i, "genbank")) 
            for warning in w:
                custom_warning_handler(warning.message, warning.category, warning.filename, warning.lineno)
        if not records:
            raise ValueError("No records found in the GenBank file.")

        # Load the YAML file
        with open(yml, 'r') as yf:  
            yaml_data = yaml.safe_load(yf)

        # Annotate records using the parsed YAML data
        annotated_records = annotate_from_yaml(yaml_data, species, codon_opt, records, pct)  

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

# update this for other feature types you want to protect
    regulatory_features = [
        "promoter", "terminator", "RBS", "regulatory",
        "rep_origin", "protein_bind", "tRNA", "oriT",
        "stem_loop", "enhancer", "5'UTR", "3'UTR",
        "polyA_signal", "LTR", "misc_binding", "tRNA",
        "rRNA", "misc_RNA"
    ]

    for record in records:
        original_features = list(record.features)  # Copy original features
        record.features = []  # Clear original list to avoid duplicates
        new_features = [] 

        for feature in original_features:
            record.features.append(feature)  # Re-add original feature

            # Add @AvoidChanges and other misc annotations based on gb record's base annotations
            if feature.type in regulatory_features:
                avoid_changes_feature = SeqFeature(
                    feature.location,
                    type="misc_feature",
                    qualifiers={"label": "@AvoidChanges"}
                )
                record.features.append(avoid_changes_feature)

            elif feature.type in ["CDS", "ORF"]:
                misc_feature = SeqFeature(
                    feature.location,
                    type="misc_feature",
                    qualifiers={"label": "@EnforceTranslation"}
                )
                record.features.append(misc_feature)

                # Handling codon options with corresponding features
                if codon_opt == 'match_codon_usage':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~match_codon_usage(species={species})"}
                    )
                    record.features.append(codon_feature)

                elif codon_opt == 'use_best_codon':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~use_best_codon(species={species})"}
                    )
                    record.features.append(codon_feature)

                elif codon_opt == 'hybrid':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~use_best_codon(species={species})"}
                    )
                    record.features.append(codon_feature)

                    enforce_changes_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~EnforceChanges(minimum_percent={pct})"}
                    )
                    record.features.append(enforce_changes_feature)

                elif codon_opt == 'harmonize_rca':
                    codon_feature = SeqFeature(
                        feature.location,
                        type="misc_feature",
                        qualifiers={"label": f"~match_codon_usage(species={species})"}
                    )
                    record.features.append(codon_feature)

                # Add comment feature last
                comment_feature = SeqFeature(
                    feature.location,
                    type="misc_feature",
                    qualifiers={"note": comment}
                )
                record.features.append(comment_feature)

        # Remove duplicates
        unique_features = []
        seen = set()
        for f in record.features:
            identifier = (
                f.type, 
                str(f.location), 
                tuple((k, tuple(v) if isinstance(v, list) else v) for k, v in f.qualifiers.items())
            )
            if identifier not in seen:
                seen.add(identifier)
                unique_features.append(f)

        record.features = unique_features  # Assign unique features back to the record

        # Marker for global constraints and optimizations
#        global_marker = SeqFeature(
#            FeatureLocation(0, len(record)),
#            type="misc_feature",
#            qualifiers={"note": "__Global constraints & optimizations start here__"})
#        record.features.append(global_marker)

        # Append additional global features that are not tied directly to existing features
        for constraint in yaml_data.get("constraints", []):
            if constraint['type'] == "EnforceTranslation":  # Skip if not applicable
                continue
            label = create_feature_label("@", constraint)
            constraint_feature = SeqFeature(
                FeatureLocation(0, len(record)),
                type="misc_feature",
                qualifiers={"label": label}
            )
            new_features.append(constraint_feature)

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

    codon_tbl = "kazusa"

    if codon_opt != "none":
        codon_data_path = "data/cleaned_coco.tsv"
        codon_data = pd.read_csv(codon_data_path, sep='\t')

        # Handling specific species mapping
        if species.lower() in ["arabidopsis", "a_thaliana"]:
            species = '3702'
            print("Using TaxID 3702 for Arabidopsis thaliana")

        # Lookup species TaxID or handle predefined species
        result = lookup_taxid(species, codon_data, "gb")

        if result[0] is None:
            print(result[1])
            sys.exit()

        # Extract taxid and comment if a valid match was found
        species, comment = result
        print(comment)

    else:
        comment = species

    annotate_genbank(comment)

