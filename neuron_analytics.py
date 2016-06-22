import catmaid_interface as ci
import networkx as nx
import scipy as sp
import numpy as np
from scipy import spatial
import scipy.sparse.csgraph as csgraph
import scipy.sparse as sparse
from collections import defaultdict


def neuron_graph(id_list, proj_opts):
    # Build a networkx neuron graph from a list of skeleton ids
    g = nx.DiGraph()

    for id in id_list:
        g.add_node(id)
        g.node[id]['name'] = ci.get_neuron_name(id, proj_opts)

    cm_edges = ci.get_connectivity_graph(id_list, proj_opts)

    for e in cm_edges:
        g.add_edge(e[0], e[1], weight=sum(e[2]))

    return g


def neuron_graph_from_annotations(annotation_list, proj_opts, anno_dict=None, append_annotations=True):
    # Build a neuron graph from a human readable list of annotation strings
    if anno_dict is None:
        anno_dict = ci.get_annotation_dict(proj_opts)

    anno_id_list = list()
    for anno in annotation_list:
        try:
            anno_id_list.append(anno_dict[anno])
        except KeyError:
            print('Not a valid key: ' + anno + ' (skipping)')

    skid_list = ci.get_ids_from_annotation(anno_id_list, proj_opts)
    g = neuron_graph(skid_list, proj_opts)

    if append_annotations:
        g = append_annotation_list(
            g, annotation_list, proj_opts, anno_dict=anno_dict)

    return g

# Given a list of annotations (as a string), add a node property to each
# skeleton containing which annotations they have


def append_annotation_list(g, annotation_list, proj_opts, anno_dict=None):
    if anno_dict is None:
        anno_dict = ci.get_annotation_dict(proj_opts)

    for anno in annotation_list:
        try:
            anno_id = [anno_dict[anno]]
            for skid in ci.get_ids_from_annotation(anno_id, proj_opts):
                if skid in g.nodes():
                    if 'annotations' in g.node[skid].keys():
                        if anno not in g.node[skid]['annotations']:
                            g.node[skid]['annotations'].append(anno)
                    else:
                        g.node[skid]['annotations'] = [anno]
        except KeyError:
            print('Not a valid key: ' + anno + ' (skipping)')
    return g


def write_node_info(g, filename, delimiter=','):
    f_nodeinfo = open(filename, 'w')
    for id, node in g.nodes_iter(data=True):
        f_nodeinfo.write(str(id))
        f_nodeinfo.write(delimiter + node['name'])
        for anno in node['annotations']:
            f_nodeinfo.write(delimiter + anno)
        f_nodeinfo.write('\n')
    f_nodeinfo.close()


class SynapseListObj:

    def __init__(self, skdata, category, proj_opts):
        if category == 'pre':
            cat_ind = 0
        elif category == 'post':
            cat_ind = 1
        else:
            raise NameError('category should be pre or post only')

        self.connids = [dat[1] for dat in skdata[1] if dat[2] == cat_ind]
        self.connid2nodeid = {dat[1]: dat[0]
                              for dat in skdata[1] if dat[2] == cat_ind}
        self.connid2loc = {dat[1]: dat[3:6]
                           for dat in skdata[1] if dat[2] == cat_ind}


class InputSynapseListObj(SynapseListObj):

    def __init__(self, skdata, proj_opts):
        SynapseListObj.__init__(self, skdata, 'post', proj_opts)

    def num(self):
        return len(self.connids)


class OutputSynapseListObj(SynapseListObj):

    def __init__(self, skdata, proj_opts):
        SynapseListObj.__init__(self, skdata, 'pre', proj_opts)
        self.num_targets = {val[0]: val[1] for val in skdata[5]}

    def num(self):
        return sum(self.num_targets.values())


class NeuronObj:

    def __init__(self, skid, proj_opts):
        self.id = skid
        skdata = ci.get_skeleton_json(self.id, proj_opts)
        self.name = skdata[4]
        self.tags = skdata[3]

        self.nodeids = [nd[0] for nd in skdata[0]]
        self.nodeloc = {nd[0]: nd[3:6] for nd in skdata[0]}
        self.node2ind = {nd[0]: i for i, nd in enumerate(skdata[0])}
        self.nodeparent = {nd[0]: nd[1] for nd in skdata[0]}
        self.radius = {nd[0]: nd[6] for nd in skdata[0]}

        temp_root = [nd[0] for nd in skdata[0] if nd[1] is None]
        self.root = temp_root[0]

        self.A = sparse.dok_matrix(
            (len(self.nodeloc), len(self.nodeloc)), dtype=np.float32)
        self.Ab = sparse.dok_matrix(
            (len(self.nodeloc), len(self.nodeloc)), dtype=np.float32)
        for key in self.nodeparent.keys():
            if self.nodeparent[key] is not None:
                self.A[self.node2ind[key], self.node2ind[self.nodeparent[key]]] = spatial.distance.euclidean(
                    self.nodeloc[key], self.nodeloc[self.nodeparent[key]])
                self.Ab[self.node2ind[key], self.node2ind[
                    self.nodeparent[key]]] = 1

        self.inputs = InputSynapseListObj(skdata, proj_opts)
        self.outputs = OutputSynapseListObj(skdata, proj_opts)

    def cable_length(self):     # Length in nm, unsmoothed
        return self.A.sum()


class SynapseDict:

    def __init__(self, neuron_list):
        self.pre_skid = dict()
        self.post_skid = defaultdict(list)
        for neuron in neuron_list:
            rel_skid = neuron.id
            for connid in neuron.inputs.connids:
                self.post_skid[connid].append(rel_skid)
            for connid in neuron.outputs.connids:
                self.pre_skid[connid] = rel_skid

    def update_synapse_dict(self, neuron_list):
        for neuron in neuron_list:
            rel_skid = neuron.id
            for connid in neuron.inputs.connids:
                self.post_skid[connid].append(rel_skid)
            for connid in neuron.outputs.connids:
                self.pre_skid[connid] = rel_skid


def dist_to_root(nrn):
    # Returns distance to root for each node in nrn as an array
    D = csgraph.shortest_path(
        nrn.A.transpose(), directed=True, unweighted=False, method='D')
    return D[nrn.node2ind[nrn.root]]


def find_end_nodes(nrn):
    # Returns a list of node ids that are end nodes (have no children)
    y = np.where(nrn.Ab.sum(0) == 0)[1]
    return [nrn.nodeids[ind] for ind in y]


def minimal_paths(nrn):
    # Returns list of lists, the minimally overlapping paths from each end
    # point toward root
    D = dist_to_root(nrn)
    ids_end = find_end_nodes(nrn)

    ends_sorted = [ids_end[ind] for ind in np.argsort(
        D[[nrn.node2ind[id] for id in ids_end]])[::-1]]
    not_visited = [True] * len(nrn.nodeids)
    min_paths = []

    for start_nd in ends_sorted:
        nd = start_nd
        min_paths.append([nd])   # Start a new list with this end as a seed
        while not_visited[nrn.node2ind[nd]] and (nrn.nodeparent[nd] is not None):
            not_visited[nrn.node2ind[nd]] = False
            nd = nrn.nodeparent[nd]
            min_paths[-1].append(nd)

    return min_paths


def neuron_list(id_list, proj_opts, syns=None):
    # Given a list of ids, build a list of NeuronObj and a SynapseDict
    # associated with them.
    neurons = [NeuronObj(id, proj_opts) for id in id_list]

    if syns is None:
        syns = SynapseDict(neurons)
    else:
        syns.update_synapse_dict(neurons)

    return (neurons, syns)


def strahler_number(neuron):
    # Computes strahler number for a neuron
    paths = minimal_paths(neuron)
    sn = {}
    for nid in neuron.nodeids:
        sn[nid] = 0

    for path in paths[::-1]:
        sn[path[0]] = 1
        for ii, nid in enumerate(path[1:]):
            if sn[nid] == sn[path[ii]]:
                sn[nid] = sn[path[ii]] + 1
            else:
                sn[nid] = sn[path[ii]]

    return sn


def split_neuron_into_components(neuron, nids):
    # Return n-component list, each element is a list of node ids in the component.
    # nids is a list of child nodes that will be split from their parent node.
    Ab_sp = copy.deepcopy(neuron.Ab)

    for id in nids:
        nind = neuron.node2ind[id]
        Ab_sp[:, nind] = 0

   ncmp, cmp_label = csgraph.connected_components(Ab_sp, directed=False)

   cmps = list()
   for cmp_val in range(ncmp):
       comp_inds = np.where( cmp_label == cmp_val )
       cmps.append( [neuron.nodeids[ind] for ind in comp_inds[0]] )
