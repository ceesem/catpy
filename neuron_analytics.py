import catmaid_interface as ci
from graph_tool import *

def make_neuron_graph( id_list, proj_opts ):
    g = Graph()
    g.add_vertex( len( id_list ) )
    skids = g.new_vertex_property("int")
    names = g.new_vertex_property("string")
    num_synapses = g.new_edge_property("double")

    id2vertex = {}
    for i, id in enumerate( id_list ):
        skids[ g.vertex(i) ] = id
        names[ g.vertex(i) ] = ci.get_neuron_name( id, proj_opts )
        id2vertex[ id ] = g.vertex(i)
    g.vertex_properties['skid'] = skids
    g.vertex_properties['name'] = names

    cm_edges = get_connectivity_graph( id_list, proj_opts )

    for e in cm_edges:
        g.add_edge( id2vertex[e[0]], id2vertex[e[1]] )
        num_synapses[ g.edge(id2vertex[e[0]], id2vertex[e[1]] ) ] = sum( e[2] )
    g.edge_properties['num_synapses'] = num_synapses

    return g
