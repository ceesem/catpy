# catmaid_to_neuron: A library for moving catmaid-y objects into NEURON-ish objects
import catmaid_interface as ci
from neuron_analytics import *
import scipy as sp
import numpy as np
import scipy.sparse.csgraph as csgraph
import scipy.sparse as sparse
from neuron import h

def segment_graph_neuron( neuron ):
    # Breaks a neuron into unbranched components and describes the graph between them
    bps = find_branch_points( neuron )
    cmps, cmp_labels = split_neuron_into_components( neuron, bps, from_parent=True)

    seg_graph = {}   # Dict of cmp_lbls that gives {child: parent}
    for bp in bps:
        par_cmp = cmp_labels[bp]
        children, dat = np.nonzero( neuron.Ab[:,neuron.node2ind[bp] ] )
        for child in children:
            seg_graph[ cmp_labels[neuron.nodeids[child] ] ] = par_cmp

    return seg_graph, cmps, cmp_labels

def assign_radius_by_straher( neuron, sn2radius, update=False ):
    # Assign a radius value (in nm) based on a strahler number lookup table
    # If sn2radius dict is shorter than the actual sn,
    # all values higher get the same radius as the max sn.
    # If update is set to True, only unassigned values get changed.
    sn = strahler_number( neuron )
    max_lookup_sn = max( sn2radius )
    if update:
        for nid in neuron.radius:
            if neuron.radius[nid] < 0:
                neuron.radius[nid] = sn2radius[ min(sn[nid],max_lookup_sn) ]
    else:
        for nid in neuron.radius:
            neuron.radius[nid] = sn2radius[ min(sn[nid],max_lookup_sn) ]
    return neuron


def assign_radius_uniform( neuron, r, update=False ):
    # Assign a single radius value (in nm) to every node in the neuron
    if update:
        for nid in neuron.radius:
            if neuron.radius[nid] < 0:
                neuron.radius[nid] = r
    else:
        for nid in neuron.radius:
            neuron.radius[nid] = r
    return neuron

class SimNeuron:
    def __init__(self, neuron, define_soma=True, Ra = 100, cm = 1):
        seg_graph, cmps, cmp_labels = segment_graph_neuron( neuron )
        self.seg_graph = seg_graph
        self.cmps = cmps
        self.cmp_labels = cmp_labels
        self.create_sections( cmps, cmp_labels, neuron )
        self.connect_sections( seg_graph )
        self.add_cable_properties( Ra, cm)
        if define_soma:
            self.define_soma_at_root( neuron )

    def create_sections(self, cmps, cmp_labels, neuron):
        d = dist_to_root(neuron)
        self.sections = []
        for ii, cmp in enumerate(cmps):
            self.sections.append( h.Section( name='cmp'+str(ii), cell=self) )
            cmpinds = [neuron.node2ind[ nid ] for nid in cmp]
            ptord = np.argsort( d[ cmpinds ] )
            for pt in ptord:
                h.pt3dadd( neuron.nodeloc[ cmp[pt] ][0] * 0.001,\
                 neuron.nodeloc[ cmp[pt] ][1] * 0.001,\
                 neuron.nodeloc[ cmp[pt] ][2] * 0.001,\
                 neuron.radius[ cmp[pt] ] * 0.001,\
                 sec=self.sections[-1]
                 ) # Remember that NEURON is tied to microns, CATMAID to nm.

    def connect_sections(self, seg_graph):
        for child in seg_graph.keys():
            self.sections[child].connect( self.sections[ seg_graph[child] ](1) )

    def find_section_from_nodeid( self, nodeid ):
        # Return the section containing a given nodeid,
        # useful for adjusting specific parts of neurons
        cmp_num = self.cmp_labels[nodeid]
        return self.sections[cmp_num]

    def add_cable_properties(self, Ra, Cm):
        for section in self.sections:
            section.Ra = Ra
            section.cm = Cm

    def define_soma_at_root( self, neuron, rad = 2000, L = 2000):
        root_section = self.find_section_from_nodeid( neuron.root )
        self.sections.append( h.Section(name='soma', cell=self) )
        root_section.connect( self.sections[-1](1) )
        self.sections[-1].L = L * 0.001
        self.sections[-1].diam = 2*rad * 0.001
        h.define_shape()

    def define_subsets(self ):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.section[0])

    def add_passive_channels(self, sections, g_pas, e_pas):
        for sec in sections:
            sec.insert('pas')
            sec.g_pas = g_pas
            sec.e_pas = e_pas
