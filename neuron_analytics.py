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
def neuron_graph_from_annotations( annotation_list, proj_opts, anno_dict = None, append_annotations=True):
    if anno_dict is None:
        anno_dict = ci.get_annotation_dict( proj_opts )

    anno_id_list = list()
    for anno in annotation_list:
        try:
            anno_id_list.append( anno_dict[anno] )
        except KeyError:
            print( 'Not a valid key: ' + anno + ' (skipping)')

    skid_list = ci.get_ids_from_annotation( anno_id_list, proj_opts )
    g = neuron_graph( skid_list, proj_opts )

    if append_annotations:
        g = append_annotation_list(g, annotation_list, proj_opts, anno_dict=anno_dict)

    return g

# Given a list of annotations (as a string), add a node property to each skeleton containing which annotations they have
def append_annotation_list( g, annotation_list, proj_opts, anno_dict = None):
    if anno_dict is None:
        anno_dict = ci.get_annotation_dict( proj_opts )

    for anno in annotation_list:
        try:
            anno_id = [ anno_dict[anno] ]
            for skid in ci.get_ids_from_annotation( anno_id, proj_opts ):
                if skid in g.nodes():
                    if 'annotations' in g.node[skid].keys():
                        if anno not in g.node[skid]['annotations']:
                            g.node[skid]['annotations'].append(anno)
                    else:
                        g.node[skid]['annotations'] = [anno]
        except KeyError:
            print( 'Not a valid key: ' + anno + ' (skipping)')
    return g

def write_node_info( g, filename, delimiter=',' ):
    f_nodeinfo = open(filename,'w')
    for id, node in g.nodes_iter(data=True):
        f_nodeinfo.write( str(id) )
        f_nodeinfo.write(delimiter+node['name'])
        for anno in node['annotations']:
            f_nodeinfo.write(delimiter+anno)
        f_nodeinfo.write('\n')
    f_nodeinfo.close()


#
# class neuron_obj:
#     def __init__( skid, proj_opts ):
#         self.id = skid;
#         skdata = ci.get_skeleton_json( self.id, proj_opts, withtags=False )
