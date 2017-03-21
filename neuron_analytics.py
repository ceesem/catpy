import catmaid_interface as ci
import scipy as sp
import numpy as np
from scipy import spatial
import scipy.sparse.csgraph as csgraph
import scipy.sparse as sparse
from collections import defaultdict
import copy

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
    # Build a neuron graph from a readable list of annotation strings
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


def append_annotation_list(g, annotation_list, proj_opts, anno_dict=None):
    # Given a list of annotations (as a string), add a node property to each
    # skeleton containing which annotations they have
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
    def __init__(self, conndata ):
        self.conn_ids = [ dat[0] for dat in conndata ]

class InputSynapseListObj(SynapseListObj):
    def __init__(self, conndata, self_id):
        SynapseListObj.__init__(self, conndata )
        self.target_node_ids = {}
        for dat in conndata:
            for skid, nid in zip( dat[1]['postsynaptic_to'], dat[1]['postsynaptic_to_node'] ):
                if skid == self_id:
                    if dat[0] in self.target_node_ids:
                        self.target_node_ids[ dat[0] ].append(nid)
                    else:
                        self.target_node_ids[ dat[0] ] = [nid]
        self.from_ids = { dat[0] : dat[1]['presynaptic_to'] for dat in conndata }
        self.from_node_ids = { dat[0] : dat[1]['presynaptic_to_node'] for dat in conndata }
    def num(self):
        return sum( map( lambda x: len(x), self.target_node_ids.values() ) )

class OutputSynapseListObj(SynapseListObj):
    def __init__(self, conndata, self_id):
        SynapseListObj.__init__(self, conndata )
        self.from_node_ids = {dat[0] : dat[1]['presynaptic_to_node'] for dat in conndata }
        self.target_ids = { dat[0] : dat[1]['postsynaptic_to'] for dat in conndata }
        self.target_node_ids = {dat[0] : dat[1]['postsynaptic_to_node'] for dat in conndata }

    def num_targets( self ):
        return {id : len(self.target_ids[id]) for id in self.target_ids}

    def num_targets_connector( self, conn_id):
        if conn_id in self.target_ids:
            return len( self.target_ids[conn_id])
        else:
            print( 'No such presynaptic connector id in neuron' )

    def num(self):
        return sum( self.num_targets().values() )

class NeuronObj:

    def __init__(self, skid, proj_opts):
        self.id = skid
        skdata = ci.get_skeleton_json(self.id, proj_opts)
        self.name = skdata[4]
        self.tags = skdata[2]

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

        pre_conn_ids = [dat[1] for dat in skdata[1] if dat[2] == 0]
        post_conn_ids = [dat[1] for dat in skdata[1] if dat[2] == 1]
        self.inputs = InputSynapseListObj( ci.get_connector_data( post_conn_ids, proj_opts ), self.id )
        self.outputs = OutputSynapseListObj( ci.get_connector_data( pre_conn_ids, proj_opts ), self.id )

    def __str__( self ):
        return self.name

    def cable_length(self):     # Length in nm, unsmoothed
        return self.A.sum()

    def get_url_to_node( self, nodeid, proj_opts, printflag=True, zoomlevel = 0 ):
        ur = ci.get_catmaid_url(
                proj_opts,
                self.nodeloc[nodeid],
                nodeid=nodeid,
                skid = self.id,
                zoomlevel=zoomlevel)
        if printflag:
            print( ur )
            return
        else:
            return ur

    def find_end_nodes(self):
        # Returns a list of node ids that are end nodes (have no children)
        y = np.where(self.Ab.sum(0) == 0)[1]
        return [self.nodeids[ind] for ind in y]

    def find_branch_points(self):
        # Returns a list of node ids that are branch points (have multiple children)
        y = np.where(self.Ab.sum(0) > 1)[1]
        return [self.nodeids[ind] for ind in y]

    def minimal_paths(self):
        # Returns list of lists, the minimally overlapping paths from each end
        # point toward root
        D = dist_to_root(self)
        ids_end = self.find_end_nodes()

        ends_sorted = [ids_end[ind] for ind in np.argsort(
            D[[self.node2ind[id] for id in ids_end]])[::-1]]
        not_visited = [True] * len(self.nodeids)
        min_paths = []

        for start_nd in ends_sorted:
            nd = start_nd
            min_paths.append([nd])   # Start a new list with this end as a seed
            while not_visited[self.node2ind[nd]] and (self.nodeparent[nd] is not None):
                not_visited[self.node2ind[nd]] = False
                nd = self.nodeparent[nd]
                min_paths[-1].append(nd)

        return min_paths

    def strahler_number(self):
        # Computes strahler number for a neuron
        paths = self.minimal_paths()
        sn = {}
        for nid in self.nodeids:
            sn[nid] = 0

        for path in paths[::-1]:
            sn[path[0]] = 1
            for ii, nid in enumerate(path[1:]):
                if sn[nid] == sn[path[ii]]:
                    sn[nid] = sn[path[ii]] + 1
                else:
                    sn[nid] = sn[path[ii]]
        return sn

    def split_into_components(self, nids, from_parent=True):
        # Return n-component list, each element is a list of node ids in the component.
        # nids is a list of child nodes that will be split from their parent node.
        # if from_parent is toggled false, parents divorce childen and not the
        # default.
        Ab_sp = copy.deepcopy(self.Ab)

        if from_parent:
            for id in nids:
                nind = self.node2ind[id]
                Ab_sp[:, nind] = 0
        else:
            for id in nids:
                nind = self.node2ind[id]
                Ab_sp[nind, :] = 0

        ncmp, cmp_label = csgraph.connected_components(Ab_sp, directed=False)

        cmps = list()
        for cmp_val in range(ncmp):
            comp_inds = np.where(cmp_label == cmp_val)
            cmps.append([self.nodeids[ind] for ind in comp_inds[0]])

        cmp_label_dict = {self.nodeids[ind]:cmp for ind,cmp in enumerate(cmp_label) }

        return cmps, cmp_label_dict

    def dist_to_root(self):
        # Returns distance to root for each node in nrn as an array
        D = csgraph.shortest_path(
            self.A.transpose(), directed=True, unweighted=False, method='D')
        return D[self.node2ind[self.root]]

    def split_by_tag( self, tag_str ):
        nids = self.tags[ tag_str ]
        cmps, cmp_label = self.split_into_components( nids )
        return cmps, cmp_label

class SynapseObject:
    def __init__(self, conn_ids, proj_opts ):
        conndata = ci.get_connector_data( post_conn_ids, proj_opts )
        self.connectors = { dat[0] : dat[1] for dat in conndata }

def neurons_from_annotations( annotation_list, proj_opts ):
    anno_dict = ci.get_annotation_dict( proj_opts )
    id_list = ci.get_ids_from_annotation( [anno_dict[anno] for anno in annotation_list], proj_opts )
    neurons = neurons_from_id_list( id_list, proj_opts )
    return neurons

def neurons_from_id_list(id_list, proj_opts ):
    # Given a list of ids, build a list of NeuronObjs
    # associated with them, if one is already pre-existing
    neurons = [NeuronObj(id, proj_opts) for id in id_list]
    return neurons

def get_adjacency_matrix( neurons, input_normalized = False ):
    # Build a weighted adjacency matrix from neurons
    A = np.zeros( (len(neurons), len(neurons)) )
    skid_to_ind = { skid:ii for ii, skid in enumerate([nrn.id for nrn in neurons])}
    ind_to_skid = { ii:skid for ii, skid in enumerate([nrn.id for nrn in neurons])}
    ids = [nrn.id for nrn in neurons]
    for nrn in neurons:
        for conn_id in nrn.outputs.target_ids:
            for targ in nrn.outputs.target_ids[conn_id]:
                if targ in ids:
                    if input_normalized is True:
                        A[ skid_to_ind[ targ ], skid_to_ind[ nrn.id ]] += 1.0 / neurons[ skid_to_ind[ targ ] ].inputs.num()
                    else:
                        A[ skid_to_ind[ targ ], skid_to_ind[ nrn.id ]] += 1
    return A, skid_to_ind, ind_to_skid

def group_adjacency_matrix( neurons, syns, groups, func=np.sum ):
    # Adjacency matrix where the entries are for groups, not neurons.
    # Groups come in a list of lists of skeleton ids.
    A, skid_to_ind, ind_to_skid = get_adjacency_matrix( neurons, syns )
    Agr = np.zeros( ( len(groups), len(groups) ) )
    for ii, grp_post in enumerate( groups ):
        for jj, grp_pre in enumerate( groups ):
            Ared = A[ [ skid_to_ind[ post ] for post in grp_post],:][:,[skid_to_ind[pre] for pre in grp_pre] ]
            Agr[ ii, jj ] = func( Ared )
    return Agr

def sort_neurons_by( neurons, sort_vector ):
    if len( sort_vector ) != len( neurons ):
        print( 'Vector must be same length as neurons' )
        return -1
    new_vector, new_neurons = zip( *sorted( zip( sort_vector,neurons ) ) )
    return new_neurons

def find_ids_by_name( neurons, name_pattern ):
    # Use regex to find sk_ids of neurons that match a given search pattern.
    ids_found = []
    return [ nrn.id for nrn in neurons if re.search(name_pattern, nrn.name) is not None ]

def number_inputs( neurons ):
    return [nrn.inputs.num() for nrn in neurons]

def number_outputs( neurons ):
    return [nrn.outputs.num() for nrn in neurons]
