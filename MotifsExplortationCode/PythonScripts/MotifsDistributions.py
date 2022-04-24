import os
import shutil
import math
import sys
from typing import *
from enum import Enum
import datetime
import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
from graph_tool.all import *
import pickle
from functools import reduce
from tqdm import tqdm, tqdm_notebook
import itertools
import argparse


parser = argparse.ArgumentParser(prog='MotifDistributionAnalyzer', description='')
parser.add_argument('--job_id', type=str,
                    help='job id as string')

args = parser.parse_args()

cluster_job_id = args.job_id
print(f"Job-ID recieved: {cluster_job_id}")


def euclidean(p1: Tuple[int, ...], p2: Tuple[int, ...]) -> float:
    return math.sqrt(reduce(lambda prev, c: prev + (c[0] - c[1]) ** 2, zip(p1, p2), 0))


def find_motifs(graph: Graph, max_motif_size: int = 10, max_visited: int = 10 ** 5):
    # g_adj_mat = get_adjacency_matrix(g)
    visited = np.zeros(g.num_vertices(), dtype=bool)
    visited_num = 0
    v_it = g.vertices()
    simple_motifs_found = {}

    def get_motif_rec(initial_vertex, current_vertex, curr_track, curr_depth):
        # global visited_num, visited
        if curr_depth >= max_motif_size:
            return None
        curr_vertex_neighbors = sorted(graph.get_all_neighbors(current_vertex),
                                       key=lambda v: euclidean(graph.vertex_properties['coordinates'][v],
                                                               graph.vertex_properties['coordinates'][current_vertex]))
        if initial_vertex in curr_vertex_neighbors and curr_depth > 1:
            return curr_track

        possible_tracks = []
        for curr_vertex_neighbor in curr_vertex_neighbors:
            if curr_vertex_neighbor == initial_vertex or curr_vertex_neighbor in curr_track:
                continue
            track = get_motif_rec(initial_vertex, curr_vertex_neighbor,
                                  curr_track + [curr_vertex_neighbor],
                                  curr_depth + 1)
            if track is None:
                continue
            else:
                possible_tracks.append(track)
        if not possible_tracks:
            return None
        else:
            return sorted(possible_tracks, key=len)[0]

    # adj_row_idx = 0
    with tqdm(total=min(max_visited, graph.num_vertices()), file=sys.stdout) as pbar:
        while visited_num < max_visited:
            pbar.set_description('processed: %d' % (1 + visited_num))

            try:
                initial_vertex = graph.vertex_index[next(v_it)]
                while visited[initial_vertex]:
                    initial_vertex = graph.vertex_index[next(v_it)]
            except StopIteration:
                break
            visited[initial_vertex] = True
            initial_v_neighbors = sorted(graph.get_all_neighbors(initial_vertex),
                                         key=lambda v: euclidean(graph.vertex_properties['coordinates'][v],
                                                                 graph.vertex_properties['coordinates'][
                                                                     initial_vertex]))
            for v_initial_neighbor in initial_v_neighbors:
                if visited[v_initial_neighbor]:
                    continue
                # v_initial_neighbor_neighbors = graph.get_all_neighbors(v_initial_neighbor)
                initial_v_motif = get_motif_rec(initial_vertex, v_initial_neighbor, [v_initial_neighbor], 1)
                if initial_v_motif is not None:
                    initial_v_motif += [graph.vertex_index[initial_vertex]]
                    initial_motif_len = len(initial_v_motif)
                    # MOTIFS_FOUND[initial_motif_len] += 1
                    simple_motifs_found[initial_motif_len] = simple_motifs_found.get(initial_motif_len, []) + [
                        initial_v_motif]

                    for v_idx in initial_v_motif:
                        visited[v_idx] = True
                    visited_num += initial_motif_len
                    pbar.update(initial_motif_len)

    return simple_motifs_found

def get_motif_properties(graph: Graph, motif_vertices_indices: List[int],
                         edge_properties: Tuple[str] = ('length', 'radii'),
                         vertex_properties: Tuple[str] = ('coordinates', 'coordinates_atlas')):
    motif_edge_properties = {}
    motif_vertex_properties = {}
    for v_prop in vertex_properties:
        motif_vertex_properties[v_prop] = [graph.vertex_properties[v_prop][v]
                                           for v in motif_vertices_indices]
    for e_prop in edge_properties:
        motif_vertex_properties[e_prop] = [graph.edge_properties[e_prop][graph.edge(*e)]
                                           for e in itertools.combinations(motif_vertices_indices, 2) if graph.edge(*e) in graph.edges()]

    # vertices_atlas = [graph.vertex_properties['coordinates_atlas'][v] for v in motif_vertices_indices]
    # vertices_artery = [graph.vertex_properties['artery_binary'][v] for v in motif_vertices_indices]
    # edge_radiis = [graph.edge_properties['radii'][e] for e in itertools.combinations(motif_vertices_indices, 2)]
    # edge_lengths = [graph.edge_properties['length'][(motif_vertices_indices[v_idx], motif_vertices_indices[v_idx+1])] for v_idx in range(len(motif_vertices_indices[:-1]))]

    # todo: get all edges between motif vertices.
    # todo: define all spatial and morphological properties as a key for each motif.
    # todo: define spatial motifs definition as a tolerance based definition.
    return motif_edge_properties, motif_vertex_properties

ssd_dir_path = f'/scratch/yishaiaz@auth.ad.bgu.ac.il/{cluster_job_id}'

ssd_path_to_graph_gt = os.path.join(ssd_dir_path, 'data_graph.gt')
print('Started loading graph from SSD')
_start = datetime.datetime.now()
g = load_graph(ssd_path_to_graph_gt)
_end = datetime.datetime.now()
t = _end - _start
print(f"Finished reading the graph in {int(t.total_seconds())} seconds")
print(f"Finished reading the graph in {int(t.total_seconds()/60)} minutes")

print('Started gathering simple motifs v_indices')
max_to_visit, max_size_of_motif = float('inf'), 10
simple_motifs = find_motifs(graph=g, max_motif_size=max_size_of_motif, max_visited=max_to_visit)
pickle.dump(simple_motifs, open(os.path.join(ssd_dir_path,
                                             f"{'all_' if max_to_visit == float('inf') else f'max_{max_to_visit}_'}simple_motifs_max_size_{max_to_visit}"),
                                'wb'), protocol=5)

print('Finished gathering simple motifs v_indices')




motifs_with_properties = {}
print('Started gathering motifs properties')
for motif_len, motifs_v_indices in tqdm(simple_motifs.items()):
    print(f"gathering motifs' properties with len={motif_len}")
    motifs_with_properties[motif_len] = motifs_with_properties.get(motif_len, []) +\
                                        [get_motif_properties(g, motif_v_indices)
                                         for motif_v_indices in motifs_v_indices]

pickle.dump(motifs_with_properties, open(os.path.join(ssd_dir_path, f"{'all_' if max_to_visit ==float('inf') else f'max_{max_to_visit}_'}motifs_with_properties_max_size_{max_to_visit}"), 'wb'), protocol=5)
print('DONE')