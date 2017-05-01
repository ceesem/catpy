import catmaid_interface as ci
import scipy as sp
import numpy as np
from scipy import spatial
import scipy.sparse.csgraph as csgraph
import scipy.sparse as sparse
from collections import defaultdict
import copy
import json

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
    def __init__(self, locs, conndata ):
        self.conn_ids = [ dat[0] for dat in conndata ]
        self.locs = locs

class InputSynapseListObj(SynapseListObj):
    def __init__(self, conndata, locs, self_id):
        SynapseListObj.__init__(self, locs, conndata )
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
    def __init__(self, conndata, locs, self_id):
        SynapseListObj.__init__(self, locs, conndata )
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

class SynapseObject:
    def __init__(self, conn_ids, proj_opts ):
        conndata = ci.get_connector_data( conn_ids, proj_opts )
        self.connectors = { dat[0] : dat[1] for dat in conndata }

def annotations_from_neurons( neurons, proj_opts ):
    anno_dat = ci.get_annotations( [nrn for nrn in neurons], proj_opts )
    anno_dict = {}
    for anno_id in anno_dat['annotations']:
        anno_dict[ int(anno_id) ] = {'str' : anno_dat['annotations'][anno_id], 'skids': [] }
    for skid in anno_dat['skeletons']:
        for anno_info in anno_dat['skeletons'][skid]:
            anno_dict[ anno_info['id'] ][ 'skids' ].append( int(skid) )
    return anno_dict

def neurons_from_annotations( annotation_list, proj_opts ):
    anno_dict = ci.get_annotation_dict( proj_opts )
    id_list = ci.get_ids_from_annotation( [anno_dict[anno] for anno in annotation_list], proj_opts )
    neurons = neurons_from_id_list( id_list, proj_opts )
    return neurons

def neurons_from_id_list(id_list, proj_opts ):
    # Given a list of ids, build a list of NeuronObjs
    # associated with them, if one is already pre-existing
    neurons = { id : NeuronObj.from_catmaid(id, proj_opts) for id in id_list}
    return neurons

def get_adjacency_matrix( neurons, input_normalized = False ):
    # Build a weighted adjacency matrix from neurons
    A = np.zeros( (len(neurons), len(neurons)) )
    skid_to_ind = { skid:ii for ii, skid in enumerate(neurons) }
    ind_to_skid = { ii:skid for ii, skid in enumerate(neurons)}
    for nrn in neurons.values():
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

# def sort_neurons_by( neurons, sort_vector ):
#     if len( sort_vector ) != len( neurons ):
#         print( 'Vector must be same length as neurons' )
#         return -1
#     new_vector, new_neurons = zip( *sorted( zip( sort_vector,neurons ) ) )
#     return new_neurons

def find_ids_by_name( neurons, name_pattern ):
    # Use regex to find sk_ids of neurons that match a given search pattern.
    ids_found = []
    return [ nrn.id for nrn in neurons if re.search(name_pattern, nrn.name) is not None ]

def number_inputs( neurons ):
    return [nrn.inputs.num() for nrn in neurons]

def number_outputs( neurons ):
    return [nrn.outputs.num() for nrn in neurons]

class NeuronObj:
    def __init__(self, neuron_info_dict):
        self.id = neuron_info_dict['id']
        self.name = neuron_info_dict['name']
        self.tags = neuron_info_dict['tags']

        self.nodeids = neuron_info_dict['nodeids']
        self.nodeloc = neuron_info_dict['nodeloc']
        self.node2ind = { nid: i for i, nid in enumerate( neuron_info_dict['nodeids'] ) }
        self.nodeparent = neuron_info_dict['nodeparent']
        self.radius = neuron_info_dict['radius']

        temp_root = [nid for nid in neuron_info_dict['nodeparent'] if neuron_info_dict['nodeparent'][nid] is None]
        self.root = temp_root[0]

        self.A = sparse.dok_matrix(
            ( len( neuron_info_dict['nodeloc'] ), len( neuron_info_dict['nodeloc'] ) ), dtype=np.float32 )
        self.Ab = sparse.dok_matrix(
            ( len( neuron_info_dict['nodeloc'] ), len( neuron_info_dict['nodeloc'] ) ), dtype=np.float32 )
        for key in neuron_info_dict['nodeparent'].keys():
            if neuron_info_dict['nodeparent'][key] is not None:
                self.A[
                    self.node2ind[ key ],
                    self.node2ind[ neuron_info_dict['nodeparent'][ key ] ]
                    ] = spatial.distance.euclidean( neuron_info_dict['nodeloc'][ key ], neuron_info_dict['nodeloc'][ neuron_info_dict['nodeparent'][ key ] ] )
                self.Ab[
                    self.node2ind[ key ],
                    self.node2ind[ neuron_info_dict['nodeparent'][ key ] ]
                    ] = 1

        self.inputs = neuron_info_dict['inputs']
        self.outputs = neuron_info_dict['outputs']

    @classmethod
    def from_catmaid(cls, skid, proj_opts):
        neuron_info_dict = {}
        neuron_info_dict['id'] = skid

        skdata = ci.get_skeleton_json( skid, proj_opts)
        neuron_info_dict['name'] = skdata[4]
        neuron_info_dict['tags'] = skdata[2]

        neuron_info_dict['nodeids'] = [nd[0] for nd in skdata[0]]
        neuron_info_dict['nodeloc'] = {nd[0]: nd[3:6] for nd in skdata[0]}
        neuron_info_dict['nodeparent'] = {nd[0]: nd[1] for nd in skdata[0]}
        neuron_info_dict['radius'] = {nd[0]: nd[6] for nd in skdata[0]}

        pre_conn_ids = [dat[1] for dat in skdata[1] if dat[2] == 0]
        post_conn_ids = [dat[1] for dat in skdata[1] if dat[2] == 1]
        conn_locs = {conn_row[1]: conn_row[3:6] for conn_row in skdata[1]}

        neuron_info_dict['inputs'] = InputSynapseListObj( ci.get_connector_data( post_conn_ids, proj_opts ), conn_locs , skid )
        neuron_info_dict['outputs'] = OutputSynapseListObj( ci.get_connector_data( pre_conn_ids, proj_opts ), conn_locs, skid)
        return cls( neuron_info_dict )

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

def msc_json_import( filename ):
    with open( filename, 'r' ) as f:
        project_data = json.load( f )

    proj_name = project_data['project_name']
    export_date = project_data['export_date']

    neurons_info_list = project_data_to_neuron_info( project_data )

    neurons = { neuron_info['id']: NeuronObj( neuron_info ) for neuron_info in neurons_info_list }
    return proj_name, export_date, neurons

def project_data_to_neuron_info( project_data ):
    # Neuron data vs neuron info is super confusing, change.
    neuron_info_list = []
    for id_str in project_data['neuron_info']:
        neuron_data = project_data['neuron_info'][id_str]
        neuron_info = { 'id': int( id_str ),
                        'name': neuron_data['name'],
                        'tags': neuron_data['node_annotations']
                        }
        neuron_info['nodeids'] = []
        neuron_info['nodeloc'] = {}
        for row in neuron_data['morphology']['node_properties']:
            neuron_info['nodeids'].append(row[0])
            neuron_info['nodeloc'][ row[0] ] = row[1:]

        neuron_info['nodeparent'] = {}
        for row in neuron_data['topology']['edge_list']:
            neuron_info['nodeparent'][ row[0] ] = row[1]

        neuron_info['radius'] = {}
        if neuron_data['morphology']['node_property_list'] is not None:
            if 'radius' in neuron_data['morphology']['node_property_list']:
                rad_ind = neuron_data['morphology']['node_property_list'].index('radius')
                for row in neuron_data['morphology']['node_properties']:
                    if row[rad_ind] > 0:
                        neuron_info['radius'][ row[0] ] = row[rad_ind]
                    else:
                        neuron_info['radius'][ row[0] ] = None
            else:
                for row in neuron_data['morphology']['node_properties']:
                    neuron_info['radius'][ row[0] ] = None
        else:
            for row in neuron_data['morphology']['node_properties']:
                neuron_info['radius'][ row[0] ] = None

        conn_post_ids = [ row[1] for row in neuron_data['connectivity_post'] ]
        connector_data_post, connector_locs_post = dict_to_connector_data( conn_post_ids, project_data )

        conn_pre_ids = [ row[1] for row in neuron_data['connectivity_pre'] ]
        connector_data_pre, connector_locs_pre = dict_to_connector_data( conn_pre_ids, project_data )

        neuron_info['inputs'] = InputSynapseListObj( connector_data_post, connector_locs_post, neuron_info['id'] )
        neuron_info['outputs'] = OutputSynapseListObj( connector_data_pre, connector_locs_pre, neuron_info['id'] )

        neuron_info_list.append( neuron_info )

    return neuron_info_list

def dict_to_connector_data( connector_list, project_data ):
    connector_data = []
    connector_locs = {}

    for conn_id in connector_list:
        conn_dat = project_data[ 'connector_info' ][ str( conn_id ) ]
        connector_locs[conn_id] = conn_dat['node_properties'][1:]
        conn_dict = {}
        conn_dict['postsynaptic_to'] = [ row[1] for row in conn_dat['connectivity_post'] ]
        conn_dict['postsynaptic_to_node'] = [ row[0] for row in conn_dat['connectivity_post'] ]
        conn_dict['presynaptic_to'] = conn_dat['connectivity_pre'][1]
        conn_dict['presynaptic_to_node'] = conn_dat['connectivity_pre'][0]
        connector_data.append( [conn_id, conn_dict] )

    return connector_data, connector_locs

def msc_json_export( neurons, filename, proj_opts, project_name=None, datestr=None, morpho_columns=None, topo_columns=None ):
    if datestr is None:
        datestr = datetime.datetime.now().strftime('%c' )

    if project_name is None:
        project_name = proj_opts['baseurl']

    neuron_info = msc_neuron_info(neurons, proj_opts, morpho_columns, topo_columns)
    connector_info = msc_connector_info( neurons, proj_opts )
    annotation_info = msc_annotation_info( neurons, proj_opts )

    # Write all the info.
    f_list = open( filename, 'w')
    json.dump( {'project_name':project_name, 'export_date':datestr, 'neuron_info':neuron_info, 'connector_info':connector_info, 'annotation_info':annotation_info}, f_list  )
    f_list.close()

def msc_connector_info( neurons, proj_opts ):
    connector_list = []
    connector_location_list = {}
    for neuron in neurons.values():
        connector_list = connector_list + neuron.inputs.conn_ids + neuron.outputs.conn_ids
        connector_location_list = {**connector_location_list,**neuron.inputs.locs,**neuron.outputs.locs}
    connector_list = list(set(connector_list))

    connector_info = {}

    conn_dat = ci.get_connector_data( connector_list, proj_opts )

    for conn in conn_dat:
        connector_info[ conn[0] ] = {}
        connector_info[ conn[0] ]['connector_annotation'] = None
        connector_info[ conn[0] ]['node_properties'] = [conn[0]] + connector_location_list[conn[0]]
        connector_info[ conn[0] ]['topology'] = None
        connector_info[ conn[0] ]['connectivity_pre'] = (conn[1]['presynaptic_to_node'], conn[1]['presynaptic_to'])
        connector_info[ conn[0] ]['connectivity_post'] =  list( zip(conn[1]['postsynaptic_to_node'], conn[1]['postsynaptic_to']) )
    return connector_info

def msc_neuron_info(neurons, proj_opts, morpho_columns=None, topo_columns=None):
    neuron_info = {}
    for nid in neurons:
        neuron_info[nid] = {}                        # Unique object id for the neuron.
        neuron = neurons[nid]
        neuron_info[nid]['name'] = neuron.name                 # Neuron name, string
        nrn_annos = na.annotations_from_neurons( {nid: neurons[nid]}, proj_opts)
        neuron_info[nid]['neuron_annotations'] = list( nrn_annos.keys() )   # A list of annotation ids

        # Morphology is everything related to nodes
        if morpho_columns is not None:                    # Morpho columns needs to be dict of dicts. First, property name, second node id.
            node_properties = [ [id]+neuron.nodeloc[id]+[morpho_columns[prop][id] for prop in morpho_columns ] for id in neuron.nodeloc ]
            node_property_list = [prop for prop in morpho_columns]
        else:
            node_properties = [ [id]+neuron.nodeloc[id] for id in neuron.nodeloc ]
            node_property_list = None
        neuron_info[nid]['morphology'] = {'node_properties': node_properties, 'node_property_list': node_property_list}

        # Topology is everything related to how nodes connect within a neuron
        topology = {}
        edge_list = [ [node_id, neuron.nodeparent[node_id]] for node_id in neuron.nodeparent ]
        edge_property_list = None

        # if topo_columns is not None:
        #     edge_list = [ [ node_id, neuron.nodeparent[ node_id ] ] + [ topo_columns[prop][ ( node_id, neuron.nodeparent[node_id] ) ] for prop in topo_columns] for node_id in neuron.nodeparent ]
        #     edge_property_list = [ prop for prop in topo_columns ]
        # else:
        #     edge_list = [ [node_id, neuron.nodeparent[node_id]] for node_id in neuron.nodeparent ] ]
        #     edge_property_list = None
        neuron_info[nid]['topology'] = {'edge_list': edge_list, 'edge_property_list': edge_property_list}

        # Connectivity is how synapses and other connections between neurons relate to anaotmy
        neuron_info[nid]['connectivity_pre'] = []
        for conn_id in neuron.outputs.from_node_ids:
            neuron_info[nid]['connectivity_pre'].append( ( neuron.outputs.from_node_ids[conn_id], conn_id) )
        neuron_info[nid]['connectivity_post'] = []
        for conn_id in neuron.inputs.target_node_ids:
            for node_id in neuron.inputs.target_node_ids[conn_id]:
                neuron_info[nid]['connectivity_post'].append( (node_id, conn_id) )

        neuron_info[nid]['connectivity_undirected'] = None

        neuron_info[nid]['node_annotations'] = neuron.tags
    return neuron_info

def msc_annotation_info( neurons, proj_opts ):
    anno_dat = na.annotations_from_neurons( neurons, proj_opts)
    annotation_info = {}
    for anno in anno_dat:
        annotation_info[anno] = {'anno_string': anno_dat[anno]['str'], 'ids': anno_dat[anno]['skids'] }
    return annotation_info
