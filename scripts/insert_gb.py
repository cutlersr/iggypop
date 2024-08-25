import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# This script merges two GenBank files by inserting the sequence and features 
# from a single-record insertion GenBank file into specified positions in 
# multiple records of an input GenBank file. The insertion occurs at positions 
# labeled with "insert_here" in the input file. The script also updates the 
# locus and accession names with a specified tag.

# Expected Input Formats:
# - **Input GenBank File**: Multi-record GenBank file with one or more sequences 
#   that include a "misc_feature" labeled "insert_here".
# - **Insertion GenBank File**: Single-record GenBank file containing the 
#   sequence and features to be inserted.

# The output is a new GenBank file with the merged sequences and updated 
# annotations.

def merge_genbank_files(input_file, insertion_file, output_file, tag):
    # Read the input GenBank files (multi-record)
    input_records = list(SeqIO.parse(input_file, "genbank"))
    
    # Read the insertion GenBank file (single-record expected)
    with open(insertion_file, "r") as insertion_handle:
        insertion_record = SeqIO.read(insertion_handle, "genbank")
    
    merged_records = []
    
    for input_record in input_records:
        # Update the locus and accession with the tag
        if not tag:
            tag = "_ins"
        input_record.name += tag
        input_record.id += tag
        if "accession" in input_record.annotations:
            input_record.annotations["accession"] = [
                acc + tag for acc in input_record.annotations["accession"]
            ]
        
        # Find the insert_here position
        insert_pos = None
        for feature in input_record.features:
            if (
                feature.type == "misc_feature" and 
                feature.qualifiers.get("label") == ["insert_here"]
            ):
                insert_pos = feature.location.end
                break
        
        if insert_pos is None:
            raise ValueError(
                "No 'insert_here' label found in the input GenBank file."
            )
        
        # Break features spanning the insert_here position
        new_features = []
        for feature in input_record.features:
            if feature.location.start < insert_pos < feature.location.end:
                # Break the feature into two
                new_features.append(
                    SeqFeature(
                        FeatureLocation(feature.location.start, insert_pos),
                        type=feature.type, qualifiers=feature.qualifiers
                    )
                )
                new_qualifiers = feature.qualifiers.copy()
                new_qualifiers["note"] = "Continued after insertion"
                new_features.append(
                    SeqFeature(
                        FeatureLocation(
                            insert_pos + len(insertion_record.seq), 
                            feature.location.end + len(insertion_record.seq)
                        ),
                        type=feature.type, qualifiers=new_qualifiers
                    )
                )
            elif feature.location.end <= insert_pos:
                new_features.append(feature)
            else:
                new_features.append(
                    SeqFeature(
                        FeatureLocation(
                            feature.location.start + len(insertion_record.seq), 
                            feature.location.end + len(insertion_record.seq)
                        ),
                        type=feature.type, qualifiers=feature.qualifiers
                    )
                )
        
        # Insert the insertion sequence and its features
        new_seq = (
            input_record.seq[:insert_pos] + 
            insertion_record.seq + 
            input_record.seq[insert_pos:]
        )
        
        for feature in insertion_record.features:
            adjusted_location = FeatureLocation(
                feature.location.start + insert_pos, 
                feature.location.end + insert_pos, 
                strand=feature.location.strand
            )
            new_features.append(
                SeqFeature(
                    adjusted_location, type=feature.type, 
                    qualifiers=feature.qualifiers
                )
            )
        
        # Create a new SeqRecord
        new_record = input_record[:]
        new_record.seq = new_seq
        new_record.features = new_features
        merged_records.append(new_record)
    
    # Write the new GenBank records to the output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(merged_records, output_handle, "genbank")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge GenBank files.")
    parser.add_argument("--i", required=True, help="Input GenBank file")
    parser.add_argument(
        "--insertion", required=True, help="Insertion GenBank file"
    )
    parser.add_argument("--o", required=True, help="Output GenBank file")
    parser.add_argument(
        "--tag", default="_ins", 
        help="Tag to append to locus and accession definitions"
    )
    
    args = parser.parse_args()
    merge_genbank_files(args.i, args.insertion, args.o, args.tag)
