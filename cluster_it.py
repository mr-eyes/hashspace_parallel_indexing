import rustworkx as rx
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Pairwise graph clustering')
parser.add_argument('-c', '--cutoff', type=int, required=True,
                    help="clustering threshold (0:100)%%")
parser.add_argument('-p', '--prefix', type=str,
                    required=True, help="pairwise csv file")
args = parser.parse_args()

input_prefix = args.prefix
# pairwise_file = input_prefix + "_kSpider_pairwise.tsv"
pairwise_file = input_prefix + "_FILTERED_kSpider_pairwise.tsv"
id_to_name_file = input_prefix + "_id_to_name.tsv"
CONTAINMENT_THRESHOLD = float(args.cutoff)
output = input_prefix + f"_kSpider_clusters_{CONTAINMENT_THRESHOLD}%.tsv"

# no_edges = 2726667056 #gtdb
# no_edges = 6928010548 #nasa_no_filtering
no_edges = 5371116604 #nasa

# with open(input_prefix + ".metadata") as metadata:
#     for line in metadata:
#         line = line.strip().split(',')
#         if line[0] == "edges":
#             no_edges = int(line[1])

# loading id_to_group_name
id_to_name = {}
with open(id_to_name_file) as F:
    next(F)
    for line in F:
        line = line.strip().split('\t')
        id_to_name[int(line[0])] = line[1]


# distance_col_idx = -1 # avg ani
distance_col_idx = -2 # avg cont

graph = rx.PyGraph()
nodes_indeces = graph.add_nodes_from(list(id_to_name.keys()))

batch_size = 10000000
batch_counter = 0
edges_tuples = []

print("[i] constructing graph")
with open(pairwise_file, 'r') as pairwise_tsv:
    next(pairwise_tsv)  # skip header
    for row in pairwise_tsv:
        row = row.strip().split('\t')
        seq1 = int(row[0]) - 1
        seq2 = int(row[1]) - 1
        distance = float(row[distance_col_idx]) * 100

        # don't make graph edge
        if distance < CONTAINMENT_THRESHOLD:
            continue

        if batch_counter < batch_size:
            batch_counter += 1
            edges_tuples.append((seq1, seq2, distance))
        else:
            graph.add_edges_from(edges_tuples)
            batch_counter = 0
            edges_tuples.clear()

    else:
        if len(edges_tuples):
            graph.add_edges_from(edges_tuples)

print("clustering...")
connected_components = rx.connected_components(graph)
print(f"number of clusters: {len(connected_components)}")
print("printing results")
single_components = 0
with open(output, 'w') as CLUSTERS:
    for component in connected_components:
        # uncomment to exclude single genome clusters from exporting
        # if len(component) == 1:
        #     single_components += 1
        #     continue
        named_component = [id_to_name[node + 1] for node in component]
        CLUSTERS.write(','.join(named_component) + '\n')

# print(f"skipped clusters with single node: {single_components}")
