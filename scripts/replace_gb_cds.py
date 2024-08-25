import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# This script replaces a specified CDS in a GenBank file 
# with a sequence from a FASTA file. It allows you to update the CDS at a 
# specified index, adjust feature locations accordingly, and clean up 
# annotations in the GenBank record. The script also allows for multiple 
# replacements, saving each modified record with a unique identifier.

# Expected Input Formats:
# - **GenBank File**: The file must contain annotated sequences with CDS 
#   features.
# - **FASTA File**: This file should contain the replacement CDS sequences.
# The "--cds" flag determines which CDS in the input GenBank file is replaced.

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Replace CDS in GenBank file with sequences from a FASTA file."
    )
    parser.add_argument("--i", required=True, help="Input GenBank file.")
    parser.add_argument(
        "--cds", type=int, required=True, 
        help="CDS number to replace (1-based index)."
    )
    parser.add_argument(
        "--replace", required=True, help="FASTA file with replacement CDSs."
    )
    parser.add_argument("--o", required=True, help="Output GenBank file.")
    return parser.parse_args()

def load_genbank_file(file_path):
    return SeqIO.read(file_path, "genbank")

def load_fasta_file(file_path):
    return list(SeqIO.parse(file_path, "fasta"))

def clean_annotations(features):
    cleaned_features = []
    for feature in features:
        if feature.type == "source":
            continue
        if "note" in feature.qualifiers:
            feature.qualifiers["note"] = [
                note for note in feature.qualifiers["note"] 
                if "Geneious type:" not in note
            ]
        if "modified_by" in feature.qualifiers:
            del feature.qualifiers["modified_by"]
        if "created_by" in feature.qualifiers:
            del feature.qualifiers["created_by"]
        cleaned_features.append(feature)
    return cleaned_features

def adjust_feature_locations(
        features, cds_feature, old_cds_length, new_cds_length):
    location_shift = new_cds_length - old_cds_length
    for feature in features:
        if feature.location.end == cds_feature.location.end:
            feature.location = FeatureLocation(
                feature.location.start, 
                feature.location.start + new_cds_length
            )
        elif feature.location.start >= cds_feature.location.end:
            feature.location = FeatureLocation(
                feature.location.start + location_shift, 
                feature.location.end + location_shift
            )
        elif feature.location.end > cds_feature.location.end:
            feature.location = FeatureLocation(
                feature.location.start, 
                feature.location.end + location_shift
            )

def replace_cds(genbank_record, cds_index, new_cds):
    features = genbank_record.features
    features = clean_annotations(features)
    cds_features = [f for f in features if f.type == "CDS"]
    if cds_index < 1 or cds_index > len(cds_features):
        raise ValueError(
            f"CDS index {cds_index} is out of range. There are {len(cds_features)} "
            f"CDS features."
        )
    cds_feature = cds_features[cds_index - 1]
    old_cds_length = len(cds_feature.location)
    new_cds_length = len(new_cds.seq)
    
    # Replace the CDS sequence
    genbank_record.seq = (
        genbank_record.seq[:cds_feature.location.start] + 
        new_cds.seq + 
        genbank_record.seq[cds_feature.location.end:]
    )
    
    # Extract accession and extra information
    accession, *extra_info = new_cds.id.split("|")
    extra_info = " ".join(extra_info)
    
    # Update the CDS feature with the new sequence and retain original qualifiers
    new_cds_feature = SeqFeature(
        FeatureLocation(
            cds_feature.location.start, 
            cds_feature.location.start + new_cds_length
        ),
        type="CDS",
        qualifiers={
            **cds_feature.qualifiers, "product": accession, "label": accession
        }
    )
    features[features.index(cds_feature)] = new_cds_feature
    
    # Add a miscellaneous feature for the extra information
    if extra_info:
        misc_feature = SeqFeature(
            FeatureLocation(
                cds_feature.location.start, 
                cds_feature.location.start + new_cds_length
            ),
            type="misc_feature", qualifiers={"note": extra_info}
        )
        features.append(misc_feature)
    
    # Adjust feature locations
    adjust_feature_locations(features, cds_feature, old_cds_length, new_cds_length)
    genbank_record.features = features
    return accession

def process_record(record, new_cds=None, cds_index=None):
    if new_cds and cds_index:
        replace_cds(record, cds_index, new_cds)
    record.features = clean_annotations(record.features)
    return record

def write_genbank_file(records, output_path):
    SeqIO.write(records, output_path, "genbank")

def main():
    args = parse_arguments()
    genbank_record = load_genbank_file(args.i)
    replacement_cdss = load_fasta_file(args.replace)
    
    # Extract prefix from the LOCUS field
    prefix = genbank_record.name
    
    # Process the original record for replacements
    modified_records = [genbank_record]
    for new_cds in replacement_cdss:
        modified_record = genbank_record[:]
        new_accession = replace_cds(modified_record, args.cds, new_cds)
        new_id = f"{prefix}_swap_{new_accession}"
        modified_record.id = new_id
        modified_record.name = new_id
        modified_record.description = new_id
        modified_records.append(modified_record)
    
    # Write the final records to the output file
    write_genbank_file(modified_records, args.o)
    print(f"Written: {args.o}")

if __name__ == "__main__":
    main()
