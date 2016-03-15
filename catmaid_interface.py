# catmaid_interface: A library for interacting with CATMAID projects through python
# Started by Tom Kazimiers (1/2013)
# Adapted by Albert Cardona (1/2013)
# Continued since by Casey Schneider-Mizell with a lot of assistance from Andrew Champion

import json
import requests
from collections import defaultdict

# requests.py authentication class for a CATMAID server with both token authentication and HTTP Basic authentication
class catmaid_auth_token( requests.auth.HTTPBasicAuth ):
    def __init__(self, token, authname=None, authpass=None):
        self.token = token
        super(catmaid_auth_token, self).__init__(authname, authpass)
    def __call__(self, r):
        r.headers['X-Authorization'] = 'Token {}'.format(self.token)
        return super(catmaid_auth_token, self).__call__(r)

# Project info
def set_project_opts( baseurl, project_id, token, authname = None, authpass = None):
    proj_opts = dict()
    proj_opts['baseurl'] = baseurl
    proj_opts['project_id'] = project_id
    proj_opts['token'] = token
    proj_opts['authname'] = authname
    proj_opts['authpass'] = authpass
    return proj_opts

# get_neuron_name: Given a skeleton ID, fetch the neuron name
def get_neuron_name( skeleton_id, proj_opts ):
    url = proj_opts['baseurl'] + '/{}/skeleton/{}/neuronname'.format( proj_opts['project_id'], skeleton_id)
    d = requests.get( url, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
    return d['neuronname']

# skeleton_object: Fetch full JSON info (plus a few other useful things) for a skeleton
# sk[0] is the list of skeleton nodes
# sk[1] is the list of connectors
# sk[2] is tags
# sk[3] is the skeleton id
# sk[4] is the neuron name
# sk[5] is the number of postsynaptic targets for each presynaptic connector
def get_skeleton_json( skeleton_id, proj_opts, withtags = True):
    if withtags:
        url = proj_opts['baseurl'] + '/{}/{}/1/1/compact-skeleton'.format( proj_opts['project_id'], skeleton_id)
    else:
        url = proj_opts['baseurl'] + '/{}/{}/1/0/compact-skeleton'.format( proj_opts['project_id'], skeleton_id)
    d = requests.get( url, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
    d.append( skeleton_id )
    d.append( get_neuron_name( skeleton_id, proj_opts ) )
    presyn_weight = post_synaptic_count( [conn[1] for conn in d[1] if conn[2] == 0], proj_opts )
    if presyn_weight!=None:
        d.append([ [cid, presyn_weight[cid]] for cid in presyn_weight.keys() ])
    else:
        d.append([])
    return d

# add_annotation: Add a single annotation to a list of skeleton IDs.
def add_annotation( annotation_list, id_list, proj_opts ):
    if type(annotation_list) is not list:
        raise TypeError('annotation_list must be a list even if just one element')
    if type(id_list) is not list:
        raise TypeError('id_list must be a list even if just one element')

    url = proj_opts['baseurl'] + '/{}/annotations/add'.format( proj_opts['project_id'])
    postdata = dict()
    for i, anno in enumerate(annotation_list):
        postdata['annotations[{}]'.format(i)] = anno
    for i, id in enumerate(id_list):
        postdata['skeleton_ids[{}]'.format(i)] = id
    d = requests.post( url, data = postdata, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) )
    return d

# get_annotation_dict: Gets the dict of annotation_id to annotation string for a project
def get_annotation_dict( proj_opts ):
    url = proj_opts['baseurl'] + '/{}/annotations/'.format(proj_opts['project_id'])
    all_annotations = requests.get( url, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
    anno_dict = { item['name']:item['id'] for item in all_annotations['annotations'] }
    return anno_dict

# get_ids_from_annotation: Given an annotation id, pull all skeletons with that annotation
def get_ids_from_annotation( annotation_id_list, proj_opts ):
    if type(annotation_id_list) is not list:
        raise TypeError('annotation_id_list must be a list even if just one element')
    url = proj_opts['baseurl'] + '/{}/annotations/query-targets'.format( proj_opts['project_id'] )
    skids = []
    for anno_id in annotation_id_list:
        postdata = { 'annotated_with' : anno_id }
        d = requests.post( url, data = postdata, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
        ids_returned = [ item['skeleton_ids'][0] for item in d['entities'] ]
        skids = skids + ids_returned
    return list( set( skids ) )

#  get_ids_by_nodecount: Ids of all skeletons with more nodes than min_nodes
def get_ids_by_nodecount( min_nodes, proj_opts  ):
    url = proj_opts['baseurl'] + '/{}/skeletons/'.format( proj_opts['project_id'] )
    d = requests.get( url, data = {'nodecount_gt': min_nodes}, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) )
    return d

# get_connected_skeleton_info: get skeleton ids, nodes, and number of synapses connected upstream/downstream for a list of ids.
def get_connected_skeleton_info( id_list, proj_opts  ):
    if type(id_list) is not list:
        raise TypeError('id_list must be a list even if just one element')

    url = proj_opts['baseurl'] + '/{}/skeletons/connectivity'.format( proj_opts['project_id'] )
    postdata = {'boolean_op' : 'OR'}
    for i, id in enumerate( id_list ):
        postdata['source_skeleton_ids[{}]'.format(i)] = id
    d = requests.post( url, data = postdata, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
    connected_inds = dict()
    connected_inds['incoming'] = {}
    connected_inds['outgoing'] = {}
    for id in d['incoming'].keys():
        connected_inds['incoming'][int(id)] = d['incoming'][id]
    for id in d['outgoing'].keys():
        connected_inds['outgoing'][int(id)] = d['outgoing'][id]
    return connected_inds

# Grow a list of skeleton ids along the connectivity network
def increase_id_list( id_list, proj_opts, min_pre=0, min_post=0, hops=1):
    if type(id_list) is not list:
        raise TypeError('id_list must be a list even if just one element')

    postdata = {'n_circles': hops,
        'min_pre': min_pre,
        'min_post': min_post}
    for i, id in enumerate(id_list):
        postdata[ 'skeleton_ids[{}]'.format(i) ] = id

    url = proj_opts['baseurl'] + '/{}/graph/circlesofhell'.format(proj_opts['project_id'])
    d = requests.post( url, data = postdata, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] )).json()

    for id in id_list:
        id_list.append( id )
    id_list = list(set(id_list))
    return id_list

# post_synaptic_count: Count how many postsynaptic targets are associated with each connector in a list of connectors
def post_synaptic_count( connector_list, proj_opts ):
    url = proj_opts['baseurl'] + '/{}/connector/skeletons'.format(proj_opts['project_id'])
    opts = {}
    for ind, id in enumerate(connector_list):
        opts[ 'connector_ids[{}]'.format(ind) ] = id
    d = requests.post( url, data = opts, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] )).json()
    nps = dict()
    for conn in d:
        nps[conn[0]] = len( conn[1]['postsynaptic_to'] )
    return nps

# write_skeletons_from_list: pull JSON files (plus key details) for
def write_skeletons_from_list( id_list, proj_opts ):
    for id in id_list:
        sk = get_skeleton_json( id, proj_opts )
        f_nodes = open( 'sk_' + str(id) + '.json' ,'w')
        json.dump(sk,f_nodes)
        f_nodes.close()

def get_connectivity_graph( id_list, proj_opts ):
    url = proj_opts['baseurl'] + '/{}/skeletons/confidence-compartment-subgraph'.format( proj_opts['project_id'] )
    opts = {}
    for i, id in enumerate(id_list):
        opts['skeleton_ids[{}]'.format(i)] = id
    d = requests.post( url, data = opts, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
    return d['edges']

# Retrieve basic statistics about skeletons like number of synapses and total length
def get_skeleton_statistics( skid, proj_opts ):
    url = proj_opts['baseurl'] + '/{}/skeleton/{}/statistics'.format( proj_opts['project_id'], skid)
    d = requests.get( url, auth = catmaid_auth_token( proj_opts['token'], proj_opts['authname'], proj_opts['authpass'] ) ).json()
    return d

#######

# Parse a file to a list of skeleton ids.
def import_ids(filename):
    f = open(filename,'r')
    id_list = []
    for line in f:
        if int(line) == -1:
            f.close()
            return []
        else:
            id_list.append(int(line))
    f.close()
    return id_list
