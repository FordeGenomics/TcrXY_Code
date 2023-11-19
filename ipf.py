#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Interspaced Palindrome Finder

## todo
# add fuzzy for flank seq
# speed up searching
    # issue with new search method, maybe?
# flag for upstream of features only
# upstream/downstream not in anno (from NZ.gb)
# add option to collapse matches

parser = argparse.ArgumentParser(description='Interspaced Palindrome Finder', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
optionalArgs = parser._action_groups.pop()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument('-g', '--genbank', type=str, help="Source GenBank file")
requiredArgs.add_argument('-f', '--fasta', type=str, help="Source FASTA file")
requiredArgs.add_argument('-r', '--reverse', action='store_const', default=False, const=True, help="Search for direct reverse tags")
requiredArgs.add_argument('-c', '--complimentary', action='store_const', default=False, const=True, help="Search for reverse complimentary tags")
optionalArgs.add_argument('-o', '--output', type=str, help="Output prefix. Leave blank to skip", default="")
optionalArgs.add_argument('-a', '--annotate', type=str, help="Output annotation info CSV name. Leave blank to skip", default="")
optionalArgs.add_argument('--min-flank', type=int, help="Minimum BP count for flanking sequence", default=4)
optionalArgs.add_argument('--max-flank', type=int, help="Minimum BP count for flanking sequence", default=8)
optionalArgs.add_argument('--min-intervening', type=int, help="Minimum BP count for intervening sequence", default=0)
optionalArgs.add_argument('--max-intervening', type=int, help="Minimum BP count for intervening sequence", default=100)
optionalArgs.add_argument('--max-distance', type=int, help="Max BP distance for annotation", default=200)

optionalArgs.add_argument('-v', '--verbose', action='store_const', default=False, const=True, help="Verbose output")
optionalArgs.add_argument('-q', '--quiet', action='store_const', default=False, const=True, help="Suppress log printing to standard out")
parser._action_groups.append(optionalArgs)
args = parser.parse_args()

# Define variable
genbank_file = args.genbank
fasta_file = args.fasta
output_file = args.output
annotate_file = args.annotate
max_intervening_length = args.max_intervening
min_intervening_length = args.min_intervening
max_flank_size = args.max_flank
min_flank_size = args.min_flank
max_distance = args.max_distance
max_length = max_intervening_length + (2 * max_flank_size)
min_length = (2 * min_flank_size) + min_intervening_length

reverse = args.reverse
complimentary = args.complimentary
annotate = args.annotate
verbose = args.verbose
quiet = args.quiet

output_header = ["seq_id", "flank_type", "flank_seq", "full_seq", "full_length", "flank_start_pos", "flank_end_pos", "intervening_length"]
annotate_header = ["seq_id", "flank_type", "flank_seq", "intervening_length", "flank_start_pos", "flank_end_pos", "flank_direction", "gene", "locus_tag", "feature_distance", "feature_start", "feature_end", "feature_strand"]

if genbank_file == None and fasta_file == None:
    print("ERROR. Need one of --genbank or --fasta required. Exiting")
    exit(1)

if genbank_file != None and fasta_file != None:
    print("ERROR. Only one of --genbank or --fasta supported. Exiting")
    print(genbank_file)
    exit(1)

if fasta_file != None and annotate_file != "":
    print("ERROR. Annotation only supported when using GenBank input. Exiting")
    exit(1)

if not (reverse or complimentary):
    print("ERROR. Please specify one or both of --reverse or --complimentary. Exiting")
    exit(1)

def print_out(msg, msg_type="LOG"):
    if quiet or (msg_type == "DEBUG" and not verbose):
        return
    else:
        print(f"{msg_type}: {msg}")

def getKmers(sequence, size, end_buff=True):
    if len(sequence) <= size:
        size = len(sequence)
    if end_buff: # give small chunks at the end to min_length
        sequence += ' ' * min_length
        for x in range(0, len(sequence) - size + 1):
            yield sequence[x:x+size]
    else:
        for x in range(0, len(sequence) - size + 1):
            yield sequence[x:x+size]

def find_tags(sequence):
    results = []
    for i in range(min_length, max_length + 1):
            subseq = Seq(sequence[:i])
            for i in range(min_flank_size, max_flank_size + 1):
                flank = Seq(subseq[:i])
                reverse_flank = flank[::-1]
                reverse_complement_flank = flank.reverse_complement()
                if subseq.endswith(reverse_flank) and reverse:
                    results.append(('Reverse', flank, subseq))
                if subseq.endswith(reverse_complement_flank) and complimentary:
                    results.append(('Reverse Compliment', flank, subseq))
    return results

def get_features(all_tags):
    features = []
    gb_features = []
    print_out("Enumerating GenBank features...")
    for feature in record.features:
        if feature.type != 'gene':
            continue
        g_location = feature.location
        g_qualifiers = feature.qualifiers
        gb_feature = {}
        gb_feature['gene'] = g_qualifiers.get('gene', [''])[0]
        gb_feature['locus_tag'] = g_qualifiers.get('locus_tag', [''])[0]
        gb_feature['start'] = int(g_location.start)
        gb_feature['end'] = int(g_location.end)
        gb_feature['strand'] = int(g_location.strand)
        gb_features.append(gb_feature)
    for seq_id, flank_type, flank_seq, full_seq, full_length, flank_start_pos, flank_end_pos, intervening_length in all_tags:
        feature_distance = len(record.seq)
        feature_closest = None
        flank_direction = "upstream"
        for gb_feature in gb_features:
            # upstream - flank end to feature start if 1, flank start to feature end if -1
            gb_strand = gb_feature['strand']
            if gb_strand == 1:
                dist = gb_feature['start'] - flank_end_pos
            elif gb_strand == -1:
                dist = flank_start_pos - gb_feature['end']
            else:
                print_out(f"unexpected strand... {gb_feature}")
                continue
            if  dist <= max_distance and dist >= 0:
                if dist < feature_distance:
                    feature_distance = dist
                    feature_strand = gb_feature['strand']
                    gene = gb_feature['gene']
                    locus_tag = gb_feature['locus_tag']
                    feature_start = gb_feature['start']
                    feature_end = gb_feature['end']
        if feature_closest != None:
            features.append((seq_id, flank_type, flank_seq, intervening_length, flank_start_pos, flank_end_pos, flank_direction, gene, locus_tag, feature_distance, feature_start, feature_end, feature_strand))
        flank_direction = "internal"
        for gb_feature in gb_features:
            # internal - flank start between feature start/end if 1, flank end between feature start/end if -1
            gb_strand = gb_feature['strand']
            if gb_strand == 1:
                if flank_start_pos >= gb_feature['start'] and flank_start_pos <= gb_feature['end']:
                    feature_distance = 'N/A'
                    feature_strand = gb_feature['strand']
                    gene = gb_feature['gene']
                    locus_tag = gb_feature['locus_tag']
                    feature_start = gb_feature['start']
                    feature_end = gb_feature['end']
                    features.append((seq_id, flank_type, flank_seq, intervening_length, flank_start_pos, flank_end_pos, flank_direction, gene, locus_tag, feature_distance, feature_start, feature_end, feature_strand))
            elif gb_strand == -1:
                if flank_end_pos <= gb_feature['end'] and flank_end_pos >= gb_feature['start']:
                    feature_distance = 'N/A'
                    feature_strand = gb_feature['strand']
                    gene = gb_feature['gene']
                    locus_tag = gb_feature['locus_tag']
                    feature_start = gb_feature['start']
                    feature_end = gb_feature['end']
                    features.append((seq_id, flank_type, flank_seq, intervening_length, flank_start_pos, flank_end_pos, flank_direction, gene, locus_tag, feature_distance, feature_start, feature_end, feature_strand))
            else:
                print_out(f"unexpected strand... {gb_feature}")
                continue
        feature_distance = len(record.seq)
        feature_closest = None
        flank_direction = "downstream"
        for gb_feature in gb_features:
            # downstream - flank start to feature end if 1, flank end to feature start if -1
            gb_strand = gb_feature['strand']
            if gb_strand == 1:
                dist = flank_start_pos - gb_feature['end']
            elif gb_strand == -1:
                dist = gb_feature['start'] - flank_end_pos
            else:
                print_out(f"unexpected strand... {gb_feature}")
                continue
            if  dist <= max_distance and dist >= 0:
                if dist < feature_distance:
                    feature_distance = dist
                    gene = gb_feature['gene']
                    locus_tag = gb_feature['locus_tag']
                    feature_start = gb_feature['start']
                    feature_end = gb_feature['end']
                    feature_strand = gb_feature['strand']
        if feature_closest != None:
            features.append((seq_id, flank_type, flank_seq, intervening_length, flank_start_pos, flank_end_pos, flank_direction, gene, locus_tag, feature_distance, feature_start, feature_end, feature_strand))
    return features

# parse genbank
entries = []
if genbank_file != None:
    print_out(f"Parsing input file {genbank_file}")
    for entry in SeqIO.parse(genbank_file, "gb"):
        entries.append(entry)
    records = [entries[0]]

if fasta_file != None:
    print_out(f"Parsing input file {fasta_file}")
    for entry in SeqIO.parse(fasta_file, "fasta"):
        entries.append(entry)
    records = entries

all_tags = []
tag_counts = {}
for record in records:
    print_out(f"Processing record {record.id}")
    print_out(f"{record}", "DEBUG")
    seq_id = record.id
    bp_add = 1
    for kmer in getKmers(record.seq, max_length):
        kmer = kmer.strip()
        tags = set(find_tags(kmer))
        for flank_type, flank_seq, full_seq in tags:
            if (flank_type, flank_seq) not in tag_counts.keys():
                tag_counts[(flank_type, flank_seq)] = []
            tag_counts[(flank_type, flank_seq)].append(seq_id)
            full_length = len(full_seq)
            flank_start_pos = bp_add
            flank_end_pos = bp_add + len(full_seq) - 1
            intervening_length = len(full_seq) - (2 * len(flank_seq))
            all_tags.append((seq_id, flank_type, str(flank_seq), str(full_seq), full_length, flank_start_pos, flank_end_pos, intervening_length))
            print_out(f"Found: {(seq_id, flank_type, flank_seq, full_seq, full_length, flank_start_pos, flank_end_pos, intervening_length)}", "DEBUG")
        bp_add += 1

n_r = 0
n_rc = 0
for seq_id, flank_type, flank_seq, full_seq, full_length, flank_start_pos, flank_end_pos, intervening_length in all_tags:
    if flank_type == "Reverse":
        n_r += 1
    else:
        n_rc += 1

if reverse:
    print_out(f"Found {n_r} direct reverse flanking regions")

if complimentary:
    print_out(f"Found {n_rc} reverse complementary flanking regions")

num_multiple = 0
multiple_df = pd.DataFrame(columns=['Motif', 'Count', 'Matches'])
for key in tag_counts.keys():
    if len(tag_counts[key]) > 1:
        multiple_df.loc[len(multiple_df)] = {'Motif': key, 'Count': len(tag_counts[key]), 'Matches': tag_counts[key]}
        num_multiple += 1

if num_multiple > 1:
    print_out(f"Found {num_multiple} identical tags:")

multiple_df = multiple_df.sort_values(by=['Count'], ascending = False)

for i, r in multiple_df.iterrows():
    print_out(f"{r['Motif']}: {r['Count']} matches, {r['Matches']}", "DEBUG")

if len(multiple_df) > 0:
    if output_file != "":
        csv_file = output_file + "_multiple.csv"
        print_out(f"Saving multiple flank CSV to {csv_file}")
        multiple_df.to_csv(csv_file, index=False)

if output_file != "":
    csv_file = output_file + "_info.csv"
    print_out(f"Saving flank information CSV to {csv_file}")
    df = pd.DataFrame(all_tags, columns=output_header)
    df.to_csv(csv_file, index=False)

if annotate_file != "":
    print_out(f"Generating annotations...")
    features = get_features(all_tags)
    print_out(f"Found {len(features)} annotations from {len(all_tags)} tags")
    print_out(f"Saving annotation information CSV to {annotate_file}")
    df = pd.DataFrame(features, columns=annotate_header)
    df.to_csv(annotate_file, index=False)

print_out("Processing done. Exiting")
exit(0)