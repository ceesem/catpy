import catmaid_interface as ci
import networkx as nx

# Build a networkx neuron graph from a list of skeleton ids
def neuron_graph( id_list, proj_opts ):
    g = nx.DiGraph()

    for id in id_list:
        g.add_node( id )
        g.node[id]['name'] = ci.get_neuron_name( id, proj_opts )

    cm_edges = ci.get_connectivity_graph( id_list, proj_opts )

    for e in cm_edges:
        g.add_edge( e[0], e[1], weight=sum(e[2]) )

    return g

# Build a neuron graph from a human readable list of annotation strings
def neuron_graph_from_annotations( annotation_list, proj_opts, anno_dict = None):
    if anno_dict is None:
        anno_dict = ci.get_annotation_dict( proj_opts )
    anno_id_list = [ anno_dict[ anno ] for anno in annotation_list ]
    skid_list = ci.get_ids_from_annotation( anno_id_list, proj_opts )
    return neuron_graph( skid_list, proj_opts )

#
# class neuron_obj:
#     def __init__( skid, proj_opts ):
#         self.id = skid;
#         skdata = ci.get_skeleton_json( self.id, proj_opts, withtags=False )
