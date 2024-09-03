import argparse
import yaml
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import warnings
from initialization import *
from pop_helpers import *

def annotate_genbank(comment):
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
        original_features = list(record.features)  # Copy original features
        record.features = []  # Clear original list to avoid duplicates
        new_features = [] 

        for feature in original_features:
            record.features.append(feature)  # Re-add original feature

            # Handling for various feature types with their corresponding misc_features
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
                    qualifiers={"note": comment}  # Using 'note' as the key for the comment
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

        record.features.extend(new_features)  # Finally, extend the features list with the new features

    return records


def lookup_taxid(input_value, data):
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Data should be a pandas DataFrame.")

    bypass_list = {'e coli', 's cerevisiae', 'b subtilis', 'c elegans', 'd melanogaster', 'm musculus', 'h sapiens'}
    normalized_input = input_value.replace('_', ' ').lower()
    
    if normalized_input in bypass_list:
        return input_value, f"No Taxid lookup for {input_value}."

    if isinstance(input_value, int) or input_value.isdigit():  # Ensure numeric input is handled correctly
        input_value = int(input_value)  # Convert to integer if it's a digit string
        matches = data[data['Taxid'] == input_value]
    elif "_" in input_value:
        matches = data[data['short_name'] == input_value]
    else:
        matches = data[data['Species'].str.lower() == input_value.lower()]

    if len(matches) == 0:
        return None, "No matches found. Please verify your input."
    elif len(matches) > 1 and "_" in input_value:
        species_list = matches['Species'].tolist()
        taxid_list = matches['Taxid'].tolist()
        return None, f"Multiple matches found: {list(zip(species_list, taxid_list))}. Please specify a unique species name or Taxid."
    
    taxid = matches['Taxid'].iloc[0]
    species = matches['Species'].iloc[0]
    short_name = matches['short_name'].iloc[0]
    comment = f"Taxid {taxid}: {species} ({short_name})"
    return taxid, comment

if __name__ == "__main__":
    tag, ofile, log_file_path, updated_defaults = initialize("format")
    globals().update(vars(updated_defaults))

    codon_data_path = "data/cleaned_coco.tsv"
    codon_data = pd.read_csv(codon_data_path, sep='\t')

    # Handling specific species mapping
    if species.lower() in ["arabidopsis", "a_thaliana"]:
        species = '3702'
        print("Using TaxID 3702 for Arabidopsis thaliana")

    # Lookup species TaxID or handle predefined species
    result = lookup_taxid(species, codon_data)

    if result[0] is None:
        print(result[1])
        sys.exit()

    # Extract taxid and comment if a valid match was found
    species, comment = result
    print(comment)

    # Check for specific codon table conditions
    if codon_tbl == "cocoputs": 
        print("Note: you can't use cocoputs data in gb mode; rerun with `--codon_tbl kazusa`.")
        sys.exit()

    # Proceed with annotation process if everything is valid
    annotate_genbank(comment)

