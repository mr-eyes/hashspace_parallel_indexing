from sourmash.distance_utils import containment_to_distance as to_ani
from tqdm import tqdm
import sys

if len(sys.argv) < 2:
    sys.exit("run: python extend_pairwise.py <index_prefix>")

input_prefix = sys.argv[1]

kmer_count_tsv = input_prefix + "_groupID_to_kmerCount.tsv"
pairwise_tsv = input_prefix + "_kSpider_pairwise.tsv"
new_pairwise_tsv = input_prefix + "_EXTENDED_kSpider_pairwise.tsv"
metadata_csv = input_prefix + ".metadata"

metadata_dict = {}

kSize = 31

# Read metadata
with open(metadata_csv) as CSV:
    for line in CSV:
        line = line.strip().split(',')
        metadata_dict[line[0]] = int(line[1])


print(f"processing with: scale={metadata_dict['scale']}, kSize={kSize}, #pairwise_comparisons={metadata_dict['edges']}")

id_to_kmer_count = {}
with open(kmer_count_tsv) as IN:
    for line in IN:
        line = line.strip().split('\t')
        id_to_kmer_count[int(line[0])] = int(line[1])
        

# TEMPORARY
id_to_name = {}
with open(input_prefix + "_id_to_name.tsv") as TSV:
    next(TSV)
    for line in TSV:
        line = line.strip().split()
        id_to_name[int(line[0])] = line[1]
    

with open(pairwise_tsv) as ORIGINAL, open(new_pairwise_tsv, 'w') as NEW:
    NEW.write(f"{next(ORIGINAL).strip()}\tavg_containment\tani\n")
    for or_line in tqdm(ORIGINAL, total=metadata_dict["edges"]):
        line = or_line.strip().split('\t')
        id_1 = int(line[0])
        id_2 = int(line[1])
        shared_kmers = int(line[2])
        max_containment = float(line[3])
        min_containment = shared_kmers / max(id_to_kmer_count[id_1],id_to_kmer_count[id_2])
        avg_containment = (max_containment + min_containment) / 2
        n_unique_kmers = (id_to_kmer_count[id_1] + id_to_kmer_count[id_2]) / 2
        ani = to_ani(avg_containment, kSize, metadata_dict["scale"], n_unique_kmers=n_unique_kmers).ani
        new_line = '\t'.join([id_to_name[id_1], id_to_name[id_2], shared_kmers, max_containment, avg_containment, ani])
        NEW.write(new_line + '\n')
        # NEW.write(f"{or_line.strip()}\t{avg_containment}\t{ani}\n")
        
        


