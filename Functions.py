import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
import itertools
from matplotlib import cm
from scipy.spatial import distance_matrix
import networkx as nx
import markov_clustering as mc
from numpyencoder import NumpyEncoder
import statistics
import os
import json
import inspect

def flatten(l):
    result = []
    for sublist in l:
        if isinstance(sublist, list):
            result += sublist
        else:
            result += [sublist]

    return result

def flatten_np(l):
    return np.asarray(np.concatenate(l))

def get_coordinate(x):
    xs = []
    ys = []
    zs = []
    for line in x:
        xs += [float(line[8])]
        ys += [float(line[9])]
        zs += [float(line[10])]
    
    return xs, ys, zs

def process_pdb(list_format, atom_types = 'C3', models = True, get_res = False):
    coor_atoms_C = []
    chains = []
    res_num = []
    result = []
    res = []
    l = [(0,6),(6,11),(12,16),(16,17),(17,20),(21,22),(22,26),
         (26,27),(30,37),(38,46),(46,54),(54,60),(60,66),(72,76),
          (76,78),(78,80)]
    model = ''
    num_model = 0
    for line in list_format:
        if 'MODEL' in line[:5]:
            num_model += 1
            if models == False and num_model > 1:
                break
            model = line.replace(' ','')


        if "ATOM" in line[:6].replace(" ","") and (len(line[17:20].replace(" ","")) == 1 or line[17:20].replace(" ","")[0] == "D") and atom_types in line[12:16]:
            new_line = [line[v[0]:v[1]].replace(" ","") for v in l ] + [model]
            #print(new_line)
            
            if new_line[5] + '_' + model not in chains:
                chains += [new_line[5] + '_' + model]
            coor_atoms_C += [new_line]

    if bool(chains) == 0:
        print("No chain found!")
        return False
    
    for chain in chains:
        sub_coor_atoms_C = [new_line for new_line in coor_atoms_C if new_line[5] == chain.split('_')[0]
                            and new_line[16] == ''.join(chain.split('_')[1:])]
        #print(sub_coor_atoms_C)
        result += [get_coordinate(sub_coor_atoms_C)]
        res_num += [[int(new_line[6]) for new_line in coor_atoms_C if new_line[5] == chain.split('_')[0]
                            and new_line[16] == ''.join(chain.split('_')[1:])]]

        res += [[new_line[4] for new_line in coor_atoms_C if new_line[5] == chain.split('_')[0] 
                 and new_line[16] == ''.join(chain.split('_')[1:])]]
    if get_res:
        return result, chains, res_num, res
    else:
        return result, chains, res_num

def list_to_range(l):
    l = list(set(l))
    l2 = []
    s = l[0]
    for p,v in enumerate(l):
        if p >= 1:
            if v == l[p-1] + 1:
                if p == len(l) - 1:
                    l2 += [range(s,v+1)]
                    
                continue
                
            e = l[p-1] + 1
            l2 += [range(s,e)]
            s = v
            
        if p == len(l) - 1:
            l2 += [range(s,v+1)]
    
    l2 = list(set(l2))

    l2 = sorted(l2, key=lambda x: x[0])
    
    return l2


'''def pymol_process(pred, res_num, name=None, color=None):
    if color is None:
                color = ['red', 'green', 'yellow', 'orange', 'blue', 'pink', 'cyan', 'purple', 'white', 'grey', 
                         'brown','lightblue', 'lightorange', 'lightpink', 'gold']

    label_set = list(set(pred))
    
    # Repeat the color list if the number of clusters exceeds the number of colors
    color = list(itertools.islice(itertools.cycle(color), len(label_set)))

    cmd = []
    for num, label in enumerate(label_set):
        label1 = [res_num[p] for p, v in enumerate(pred) if v == label]
        clust_name = name + f'_cluster_{num}' if name is not None else f'cluster_{num}'
        cmd.append(command_pymol(label1, clust_name, color[num]))

    return cmd'''

def generate_colors(num_colors):
    colormap = cm.get_cmap('hsv', num_colors)
    return [colormap(i) for i in range(num_colors)]

def pymol_process(pred, res_num, name=None, color=None):
    if color is None:
        color = ['red', 'green', 'yellow', 'orange', 'blue', 'pink', 'cyan', 'purple', 'white', 'grey', 
                    'brown','lightblue', 'lightorange', 'lightpink', 'gold']

    label_set = list(set(pred))

    if len(label_set) > len(color):
        # Generate additional colors dynamically
        additional_colors = generate_colors(len(label_set) - len(color))
        # Convert additional colors from RGBA to hex format
        additional_colors_names = ['{:02d}'.format(i) for i in range(len(color), len(label_set))]
        color.extend(additional_colors_names)

    cmd = []
    for num, label in enumerate(label_set):
        label1 = [res_num[p] for p, v in enumerate(pred) if v == label]
        clust_name = name + f'cluster_{num+1}' if name is not None else f'cluster_{num+1}'
        cmd.append(command_pymol(label1, clust_name, color[num]))

    return cmd

def command_pymol(l, name, color):
    l2 = list_to_range(l)
    mess = f'select {name}, res '
    for p,r in enumerate(l2):
        if len(r) > 1:
            mess += f'{r[0]}-{r[-1]}'
            if p != len(l2) - 1:
                mess += '+'
        else:
            mess += f'{r[0]}' + '+'
    mess += f'; color {color}, {name}'
    print(mess)
    return mess

def distance_2arrays(arr1, arr2):
    dist = 1
    for i in range(len(arr1)):
        if arr1[i] == arr2[i]:
            dist -= 1/len(arr1)
    
    return dist

def join_clusters(list_cluster):
    prev_result = list_cluster
    result = []
    cont = True
    while cont:
        cont = False
        for cluster1, cluster2 in itertools.combinations(prev_result, 2):
            if any(i in cluster1 for i in cluster2):
                if cluster1 in result:
                    result.remove(cluster1)
                if cluster2 in result:
                    result.remove(cluster2)
                
                result += [tuple(set(cluster1 + cluster2))]
                cont = True
        
        result = list(set(result))
        prev_result = result

    return prev_result

def merge_clusters(clusters, G, k):
    # Calculate modularity contribution for each cluster
    def modularity_gain(cluster):
        internal_edges = G.subgraph(cluster).size(weight='weight')
        total_degree = sum(dict(G.degree(cluster, weight='weight')).values())
        return internal_edges - (total_degree ** 2) / (2 * G.size(weight='weight'))

    while len(clusters) > k:
        # Find the pair of clusters whose merge gives the highest modularity gain
        best_pair = None
        best_modularity = float('-inf')
        for (i, cluster1), (j, cluster2) in itertools.combinations(enumerate(clusters), 2):
            combined_cluster = cluster1.union(cluster2)
            modularity = modularity_gain(combined_cluster)
            if modularity > best_modularity:
                best_modularity = modularity
                best_pair = (i, j)

        # Merge the best pair of clusters
        cluster1, cluster2 = best_pair
        clusters[cluster1] = clusters[cluster1].union(clusters[cluster2])
        clusters.pop(cluster2)

    return clusters

def build_graph(data, distance = 8, weight = False):
    # Compute the pairwise distance matrix
    #dist_matrix = distance_matrix(data1['data'], data1['data'])
    dist_matrix = distance_matrix(data, data)
    vectorize_normalize = np.vectorize(contact_prob)
    mod_vectorize_normalize = np.vectorize(mod_contact_prob)

    weight_matrix =  vectorize_normalize(dist_matrix)
    mod_weight_matrix = mod_vectorize_normalize(dist_matrix)
    # Create a graph
    G = nx.Graph()

    # Add nodes
    #for i, point in enumerate(data1['data']):
    for i, point in enumerate(data):
        G.add_node(i, pos=(point[0], point[1]))

    # Add edges based on the distance threshold
    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            if dist_matrix[i, j] < distance:
                if weight in ("True", "T", "true"):
                    #if j == i+1:
                    if j == i+1:
                        G.add_edge(i, j, weight = mod_weight_matrix[i,j])
                    else:
                        G.add_edge(i, j, weight = weight_matrix[i,j])
                else:
                    G.add_edge(i,j)

    return G

def cluster_algo(*args):
    data = args[0]
    G = build_graph(data, args[2], args[3])
    if args[1] == 'H':
        print("Executing Hierachical-based clustering...")
        _, result1 = top_down_commu2(G, args[4], args[6])

        result1 = [list(range(i[0], i[1])) for i in result1]

        #bot_up_commu2(G, result1)

        result = bot_up_commu2(G, result1, args[5])
        
    elif args[1] == 'M':
        print("Executing Markov clustering...")
        dist_matrix = distance_matrix(data, data)
        adj_matrix = np.array([[1 if dist_matrix[j][i] < args[2] and i != j else 0 for i in range(len(dist_matrix[j]))] 
                       for j in range(len(dist_matrix))])
        
        model = mc.run_mcl(adj_matrix, expansion=args[4], inflation=args[5])
        
        result = mc.get_clusters(model)

    elif args[1] == 'L':
        print("Executing Louvain clustering...")

        result = nx.community.louvain_communities(G, resolution =  args[4])

        #result = merge_clusters(result1, G, args[4])
        #result = merge_clusters(result1, G, 2)
    
    elif args[1] == 'C':
        print("Executing Clauset-Newman-Moore clustering...")
        result = nx.community.greedy_modularity_communities(G, resolution = args[4])

    elif args[1] == 'G':
        print("Executing Girvan-Newman clustering...")
        result1 = nx.community.girvan_newman(G)
        for _ in range(0):
            next(result1)

        result = tuple(sorted(c) for c in next(result1))
        
    else:
        print(args[1])
        sys.exit("Non recognized algorithm!")

    return result

def check_C(result, threshold):
    data = []
    if result == False or len(result) == 0:
        return False

    else:
        for t in range(len(result[0])):
            if len(result[0][t][0]) < threshold:
                continue
        
            l = [[result[0][t][0][i], result[0][t][1][i], result[0][t][2][i]] for i in range(len(result[0][t][0]))]
            data += [np.array(l)]

        return data, [i for i in result[-1] if len(i) >= threshold]

    
def domain_overlap_matrix(lists_label, list_residue = None): #Order in lists_label: ground_truth, prediction 
    if list_residue == None:
        list_residue = range(len(lists_label[0]))
    
    group_label = {'pred': lists_label[1], 'true': lists_label[0]}

    group_residue = {'pred': {}, 'true': {}}

    for key in group_label.keys():
        for label in set(group_label[key]):
            #group_residue[key][label] = list_to_range([lists_residue[key][i] for i in range(len(lists_residue[key])) if lists_label[key][i] == label])
            group_residue[key][label] = [list_residue[i] for i in range(len(list_residue)) if group_label[key][i] == label]

    domain_matrix = []
    for label in sorted(set(group_label['pred'])): 
        row = []
        for label2 in sorted(set(group_label['true'])):
            row += [len([intersect for intersect in group_residue['pred'][label] if intersect in group_residue['true'][label2]])]

        domain_matrix += [row]


    min_labels = [min(set(group_label['pred'])), min(set(group_label['true']))]
    #return domain_matrix
    return domain_matrix, min_labels

def NDO(domain_matrix, len_rnas, min_labels = [0,0]):
    #print(domain_matrix, [min_labels[0]**2, min_labels[1]**2])
    domain_matrix_no_linker = [row[min_labels[1]**2:] for row in domain_matrix[min_labels[0]**2:]]
    
    domain_matrix_no_linker = np.asarray(domain_matrix_no_linker)
    domain_matrix = np.asarray(domain_matrix)
    
    sum_col = np.sum(domain_matrix, axis = 0)
    
    sum_row =  np.sum(domain_matrix, axis = 1)
    
    max_col = np.amax(domain_matrix_no_linker, axis = 0)
    
    max_row = np.amax(domain_matrix_no_linker, axis = 1)
    
    Y = 0
    print(domain_matrix,domain_matrix_no_linker, sum_col, sum_row, max_col, max_row, sep='\n')
    for row in range(domain_matrix_no_linker.shape[0]):

        Y += 2*max_row[row] - sum_row[row+min_labels[0]**2]
        
    for col in range(domain_matrix_no_linker.shape[1]):
        Y += 2*max_col[col] - sum_col[col+min_labels[1]**2]

    score = Y*100/(2*(len_rnas - sum_col[0]*(min_labels[1])**2))
    #score = Y*100/(2*(len_rnas))
    
    return score

def domain_boundary(lists_label, list_residue = None): #Order in lists_label: ground_truth, prediction 
        if list_residue == None:
            list_residue = range(len(lists_label[0]))
        
        group_label = {'pred': lists_label[1], 'true': lists_label[0]}

        group_residue = {'pred': {}, 'true': {}}
        
        group_boundary = {'pred': {}, 'true': {}}
        
        for key in group_label.keys():
            list_boundary = []
            for label in set(group_label[key]):
                #group_residue[key][label] = list_to_range([lists_residue[key][i] for i in range(len(lists_residue[key])) if lists_label[key][i] == label])
                group_residue[key][label] = [ list_to_range(list_residue[i] for i in range(len(list_residue)) if group_label[key][i] == label)]
            
            
                list_boundary += flatten([[(j[0],j[-1]) for j in i] for i in group_residue[key][label]])
            
            group_boundary[key] = list_boundary
            
        return group_boundary

def merge_ele_matrix(mtx, list_range_true, list_range_pred):
    merged_mtx_row = []
    for row in mtx:
        new_row = row.copy()  # Create a copy of the row to avoid modifying the original row
        for label_pos in range(len(list_range_true)):
            s = label_pos
            e = s + len(list_range_true[label_pos])
            new_row = new_row[:s] + [sum(new_row[s:e])] + new_row[e:]
        merged_mtx_row.append(new_row)
    
    merged_mtx_row = np.array(merged_mtx_row).T.tolist()
    
    merged_mtx = []
    for row in merged_mtx_row:
        new_row = row.copy()  # Create a copy of the row to avoid modifying the original row
        for label_pos in range(len(list_range_pred)):
            s = label_pos
            e = s + len(list_range_pred[label_pos])
            new_row = new_row[:s] + [sum(new_row[s:e])] + new_row[e:]
        merged_mtx.append(new_row)

        
    merged_mtx = np.array(merged_mtx).T  # Transpose to get the correct shape
    
    return merged_mtx

def domain_distance(segment1, segment2):
    d = abs(min(segment1) - min(segment2))
    d += abs(max(segment1) - max(segment2))
    
    return d/2

def domain_distance_matrix(lists_label, list_residue = None):#Order in lists_label: ground_truth, prediction 
    if list_residue == None:
        list_residue = range(len(lists_label[0]))
    
    group_label = {'pred': lists_label[1], 'true': lists_label[0]}

    group_residue = {'pred': {}, 'true': {}}
    for key in group_label.keys():
        for label in set(group_label[key]):
            group_residue[key][label] = [list_residue[i] for i in range(len(list_residue)) if group_label[key][i] == label]
    
    domain_distance_mtx = []
    if len(set(group_label['pred'])) > 0:
        for label1 in sorted(set(group_label['pred'])):
            if bool(group_residue['pred'][label1]):
                for segment1 in list_to_range(group_residue['pred'][label1]):
                    row = []
                    for label2 in sorted(set(group_label['true'])):
                        for segment2 in list_to_range(group_residue['true'][label2]):
                            #print(segment1, segment2)
                            x = domain_distance(segment2, segment1)
                            row += [x]
                            #print(x, end = '\n\n')

                    domain_distance_mtx += [row]

        lst_to_range_true = [list_to_range(group_residue['true'][label1]) for label1 in sorted(set(group_label['true'])) if bool(group_residue['true'][label1])]
        lst_to_range_pred = [list_to_range(group_residue['pred'][label1]) for label1 in sorted(set(group_label['pred'])) if bool(group_residue['pred'][label1])]

    return domain_distance_mtx, lst_to_range_true, lst_to_range_pred

def DBD(domain_distance_mtx, list_range_true, list_range_pred, threshold = 50):
    scoring_mtx = []
    for row in domain_distance_mtx:
        scoring_mtx += [[threshold - i if i < threshold else 0 for i in row]]
    
    merged_scoring_mtx = merge_ele_matrix(scoring_mtx, list_range_true, list_range_pred)

    merged_scoring_mtx = np.asarray(merged_scoring_mtx)
    scoring_mtx = np.asarray(scoring_mtx)
    print(merged_scoring_mtx, scoring_mtx, sep = '\n')
    if scoring_mtx.shape[0] >= scoring_mtx.shape[1]:
        max_row = np.amax(merged_scoring_mtx, axis = 1)
        total_score = sum(max_row)/(threshold*scoring_mtx.shape[0])
    else:
        max_col = np.amax(merged_scoring_mtx, axis = 0)
        total_score = sum(max_col)
        total_score = sum(max_col)/(threshold*scoring_mtx.shape[1])
    
    return total_score

def inf(tup_of_tups, val):
    start_eles = [min(tup) for tup in tup_of_tups if min(tup) <= val]
    end_eles = [max(tup) for tup in tup_of_tups if max(tup) <= val]
    
    inf_start = max(start_eles) if bool(start_eles) != False else -1
    inf_end = max(end_eles) if bool(end_eles) != False else -1

    return inf_start, inf_end

def sup(tup_of_tups, val):
    start_eles = [min(tup) for tup in tup_of_tups if min(tup) >= val]
    end_eles = [max(tup) for tup in tup_of_tups if max(tup) >= val]
    
    sup_start = min(start_eles) if bool(start_eles) != False else val
    sup_end = min(end_eles) if bool(end_eles) != False else val

    return sup_start, sup_end

def find_tuple_index(tuples, element, order = False):
    for i, t in enumerate(tuples):
        if element in t:
            if order == False:
                return i
            else:
                if order == "first":
                    if element == t[0]:
                        return i
                else:
                    if element == t[-1]:
                        return i
    return -1  # Return -1 if the element is not found in any tuple

def contact_prob(d, d0 = 8, sig = 1.5):
    p = 1/(1+math.e**((d - d0)/sig))
    
    return p

def mod_contact_prob(d, d0 = 8, sig = 1.5):
    p = 1/(1+0.5*math.e**((d - d0)/(sig)))
    
    return p

def check_clique(nodes, edges, not_connected = 0):
    edges = list(set([tuple(sorted(i)) for i in edges]))
    nodes = list(set(nodes))
    
    #connect_dict = {node: sum(1 for edge in edges if node in edge) >= len(nodes) - 1 - not_connected for node in nodes}
    for node in nodes:
        node_count = sum(1 for edge in edges if node in edge)
        if node_count < len(nodes) - 1 - not_connected:
            return False
    
    return True

def mod_modularity_score(nodes, edges, alpha = 1):
    edges = list(set([tuple(sorted(i)) for i in edges]))
    nodes = sorted(list(set(nodes)))
    
    degrees = [len([i for i in edges if node in i]) for node in nodes]
    print(len(nodes), len(edges))
    score = 0
    for inode1, inode2 in itertools.combinations(range(len(nodes)), 2):
        #A = 1 if (nodes[inode1],nodes[inode2]) in edges else 0
        
        #print(A)
        #score += A - degrees[inode1]*degrees[inode2]/len(edges)
        x = degrees[inode1]*degrees[inode2]/(2*len(edges)*(len(nodes)-1)*alpha) 
        score += x - 1
        #print(x)
    
    #score -= 1
    
    return score

def modularity_score(nodes, edges, alpha = 1):
    edges = list(set([tuple(sorted(i)) for i in edges]))
    nodes = sorted(list(set(nodes)))
    
    degrees = [len([i for i in edges if node in i]) for node in nodes]
    score = 0
    for inode1, inode2 in itertools.combinations(range(len(nodes)), 2):
        A = 1 if (nodes[inode1],nodes[inode2]) in edges else 0
        
        score += A - degrees[inode1]*degrees[inode2]/(2*len(edges))

    score /= (2*len(edges))
    
    return score

def largest_smaller_than(lst, value):
    # Initialize variables to store the largest element found and its index
    largest = None
    largest_index = -1
    
    # Iterate through the list with index
    for index, elem in enumerate(lst):
        # Check if the element is smaller than the given value
        if elem <= value:
            # If largest is None or current element is larger than largest found so far
            if largest is None or elem > largest:
                largest = elem
                largest_index = index
    
    return largest, largest_index,value

def sort_dict_by_value(input_dict, descending=True):
    # Sort the dictionary based on values
    sorted_items = sorted(input_dict.items(), key=lambda x: x[1], reverse=descending)
    
    # Create a new dictionary with sorted items
    sorted_dict = {k: v for k, v in sorted_items}
    
    return sorted_dict

def tup_pos_process(tup_of_tup):
    result = []
    for tup in tup_of_tup:
        s = []
        for i in tup:
            s += list(range(i[0], i[1]))
        
        s = tuple(s)
        result += [s]
    
    return result

def top_down_commu(graph, min_value = 0.4):
    nodes = sorted(list(graph.nodes).copy())
    edges = list(graph.edges).copy()
    
    # Each subgraph is corresponse to each segment 
    subgraphs = [graph.copy()]
    list_pos = [(0,max(nodes)+1)]
    p = 0
    while True:
        s = 0
        temp_list_pos = list_pos.copy()
        temp_subgraphs = subgraphs.copy()
        for pos, segment in enumerate(list_pos):
            dict_score = {}
            nodes = list(temp_subgraphs[pos+s].nodes)
            for node1 in nodes[min(segment) + 30: max(segment)-30]:
                node_list1 = nodes[:nodes.index(node1)]; node_list2 = nodes[nodes.index(node1):]
                #print(node1)
                dict_score[node1] =  nx.community.modularity(subgraphs[pos], [set(node_list1), set(node_list2)])
            
            if bool(dict_score.keys()) == False:
                continue
            
            key = max(dict_score, key=dict_score.get)
            max_value = dict_score[key]
            
            if max_value < min_value:
                continue
                
            temp_list_pos = temp_list_pos[:pos+s] + [(min(segment), key), (key, max(segment))] + temp_list_pos[pos+1+s:]
            
            subgraph1 = temp_subgraphs[pos].copy(); subgraph2 = temp_subgraphs[pos].copy()
            #print(key, segment, max_value)
            subgraph1.remove_nodes_from(list(range(key, max(segment))))
            subgraph2.remove_nodes_from(list(range(min(segment), key)))            
            subgraph1.remove_edges_from([edge for edge in subgraph1.edges if any(e for e in edge if e not in list(subgraph1.nodes))])
            subgraph2.remove_edges_from([edge for edge in subgraph2.edges if any(e for e in edge if e not in list(subgraph2.nodes))])
            
            temp_subgraphs = temp_subgraphs[:pos+s] + [subgraph1, subgraph2] + temp_subgraphs[pos+1+s:]
            
            s += 1

        list_pos = temp_list_pos
        subgraphs = temp_subgraphs
        
        p += 1
        if p == 40:
            break
    
    print(list_pos)
    list_pos = [list(range(i[0], i[1])) for i in list_pos]
    return list_pos 

def bot_up_commu(graph, segment_list):
    score_list = []
    segments = sorted(segment_list.copy(), key=len)
    #segments = [list(range(i[0], i[1])) for i in segment_list]
    print(len(segments))
    #new_segments = segments.copy()
    #for segment in segments:
    p = 0
    while p <= len(segments) - 1:
        segment = segments[p]
        other_segments = [i for i in segments if i != segment]
        others = [j for i in segments for j in i if i != segment]
        
        score = nx.community.modularity(graph, [set(segment), set(others)])
        score_list += [score]    
        #print('all: ', score)
        
        max_score = -1
        max_set = []
        max_other = []
        for pos, other in enumerate(other_segments):
            set2 = set(other + segment)
            set3 = set([j for i in other_segments for j in i if i != other])
            score2 = nx.community.modularity(graph, [set2, set3])
            if score2 > max_score:
                max_score = score2
                max_set = set2
                max_other = other
                #max_pos = pos
            
        if max_score >= abs(1*score):
            #min_pos = min(p, max_pos)
            #others2 = [i for i in other_segments if i != max_other]
            #segments = others2[:min_pos] + [list(max_set)] + others2[min_pos+1:]
            segments = [list(max_set)] + [i for i in other_segments if i != max_other]
            p = 0
        else:
            p += 1

        #print(score, max_score)
    
    #print(p)
    return segments 

def top_down_commu2(graph, min_value = 0.4, resolution = 1):
    nodes = sorted(list(graph.nodes).copy())
    edges = list(graph.edges).copy()
    
    # Each subgraph is corresponse to each segment 
    subgraphs = graph.copy()
    segment = (min(nodes),max(nodes))

    dict_score = {}

    for node1 in nodes:
        node_list1 = nodes[:nodes.index(node1)]; node_list2 = nodes[nodes.index(node1):]
        #print(node1)
        dict_score[node1] =  nx.community.modularity(subgraphs, [set(node_list1), set(node_list2)], resolution = resolution)
        

    if bool(dict_score.keys()) == False:
        #print('dict_score is empty')
        return graph, [(min(nodes),max(nodes)+1)]
    
    key = max(dict_score, key=dict_score.get)
    max_value = dict_score[key]
    
    if max_value < min_value or key <= min(nodes) + 30 or key >= max(nodes) - 30:
        #print('max_value is smaller than threshold', max_value)
        return graph, [(min(nodes),max(nodes)+1)]
    
    subgraph1 = subgraphs.copy(); subgraph2 = subgraphs.copy()
    subgraph2.remove_nodes_from(list(range(min(segment), key)))  
    subgraph1.remove_nodes_from(list(range(key, max(segment)+1)))          
    subgraph1.remove_edges_from([edge for edge in subgraph1.edges if any(e for e in edge if e in list(subgraph2.nodes))])
    subgraph2.remove_edges_from([edge for edge in subgraph2.edges if any(e for e in edge if e in list(subgraph1.nodes))])

    subgraph1,list_pos1 = top_down_commu2(subgraph1, min_value)
    subgraph2,list_pos2 = top_down_commu2(subgraph2, min_value)
    
    #temp_subgraphs = temp_subgraphs[:pos+s] + [subgraph1, subgraph2] + temp_subgraphs[pos+1+s:]
    list_pos = list_pos1 + list_pos2

    return graph, list_pos 

def bot_up_commu2(graph, segment_list, ratio = 0.1):
    score_list = []
    segments = sorted(segment_list.copy(), key=len)
    #segments = [list(range(i[0], i[1])) for i in segment_list]
    #new_segments = segments.copy()
    #for segment in segments:
    p = 0
    while p <= len(segments) - 1 and len(segments) > 1:
        #print(p, len(segments))
        segment = segments[p]
        other_segments = [i for i in segments if i != segment]
        others = [j for i in segments for j in i if i != segment]
        
        max_score = 0; max_segment = []
        for p2, segment2 in enumerate(other_segments):
            subgraph = graph.copy()
            subgraph.remove_nodes_from([i for i in list(graph.nodes) if i not in segment2 + segment])
            subgraph.remove_edges_from([edge for edge in subgraph.edges if any(e for e in edge if e not in list(subgraph.nodes))])
            #score = inter_community_edges(subgraph, [segment2, segment])
            score2 = inter_community_edges2(subgraph, [segment2, segment])

            if score2 > max_score:
                max_score = score2
                max_segment = segment2
                max_subgraph = subgraph
                max_p = p2

        #threshold = 0.1/min(len(segment), len(max_segment))
        threshold = ratio*min(len(segment), len(max_segment))
        old_segments = segments.copy()  
        #print(max_score, threshold)
        if max_score >= threshold and bool(max_segment):
            segments = [list(set(max_segment + segment))] + [i for i in other_segments if i != max_segment]
            segments = sorted(segments, key=len)
            p = max(min(p, max_p),0)
            #p = 0

        else:  
            p += 1
            continue
    '''for segment1,segment2 in itertools.combinations(segments, 2):
        subgraph = graph.copy()
        subgraph.remove_nodes_from([i for i in list(graph.nodes) if i not in segment1 + segment2])
        subgraph.remove_edges_from([edge for edge in subgraph.edges if any(e for e in edge if e not in list(subgraph.nodes))])
        score2 = inter_community_edges2(subgraph, [segment1, segment2])
        threshold = 0.3*min(len(segment1), len(segment2))
        print(segments.index(segment1), segments.index(segment2),
              score2, threshold)'''
    
    return segments

def inter_community_edges(G,segments):
    betweenness = nx.edge_betweenness_centrality(G)
    scores = []
    s = 0
    for edge, score in betweenness.items():
        if ((edge[0] in segments[0])^(edge[1] in segments[0])):
            if score > 0:
                scores += [score]

    if s == 0:
        return 0
    
    normalized_scores = sum(scores)/s
    
    return normalized_scores

def inter_community_edges2(G,segments):
    nodes = []

    for edge in G.edges:
        if ((edge[0] in segments[0])^(edge[1] in segments[0])):
            if edge[0] not in nodes:
                nodes += [edge[0]]
            if edge[1] not in nodes:
                nodes += [edge[1]]
    
    return len(nodes)

def process_cluster_format(clust_lst, res_lst = None):
    if res_lst == None:
        res_lst = list(range(1,len(clust_lst)+1))

    clust_by_res = []
    set_clust = set(clust_lst)
    for clust in set_clust:
        sublst = []
        for pos, res in enumerate(res_lst):
            if clust_lst[pos] ==  clust:
                sublst += [res]

        clust_by_res += [sublst]

    return clust_by_res
    
def split_pdb_by_clusters(pdb_file, clusters, output_prefix, chain=None):
    """
    Splits a PDB file into multiple PDB files based on provided clusters of residues.
    If a chain is specified, only residues from that chain will be processed.

    Parameters:
        pdb_file (str): Path to the input PDB file.
        clusters (list of list of int): List of clusters, where each cluster is a list of residue indices.
        output_prefix (str): Prefix for output files.
        chain (str, optional): Chain ID to filter residues by. If None, all chains are processed.
    """

    # Read the PDB file
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    # Create a dictionary to store lines for each cluster
    cluster_lines = {i: [] for i in range(len(clusters))}

    # Process each line in the PDB file
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract residue sequence number and chain ID
            residue_seq = int(line[22:26].strip())  # Extract residue sequence number
            residue_chain = line[21:22].strip()  # Extract chain ID (column 22)

            # If a chain is specified, skip lines that don't match the chain
            if chain and residue_chain != chain:
                continue

            # Check which cluster this residue belongs to
            for cluster_index, cluster in enumerate(clusters):
                if residue_seq in cluster:
                    cluster_lines[cluster_index].append(line)
                    break

    # Write the output files for each cluster
    for cluster_index, cluster in cluster_lines.items():
        if cluster:
            output_file = f"{output_prefix}_cluster_{cluster_index + 1}.pdb"
            with open(output_file, 'w') as outfile:
                outfile.writelines(cluster)
            print(f"Wrote {len(cluster)} lines to {output_file}")
    