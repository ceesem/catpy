import sys
from collections import defaultdict
import os
import catmaid_interface as ci
import cardonalab_project_parser as cpp

# Given a list of annotations, get skeletons from them and perhaps their neighbors
def main():

    project_flag = sys.argv[1]
    input_args = { sys.argv[i] : sys.argv[i+1] for i in range(2,len(sys.argv),2) }

    proj_opts = cpp.project_map( project_flag )
    anno_dict = ci.make_annotation_dict(proj_opts))

    arg_ids_to_expand = [ anno_dict[arg] for arg in input_args.keys() if input_args[arg]=='1' ]
    arg_ids_noexpand = [ anno_dict[arg] for arg in input_args.keys() if input_args[arg]=='0' ]

    if len(arg_ids_to_expand)>0:
        id_list_base = ci.get_ids_from_annotation(arg_ids_to_expand, proj_opts )
    else:
        id_list_base = []

    if len(arg_ids_noexpand)>0:
        id_list_noexpand = ci.get_ids_from_annotation( arg_ids_noexpand, proj_opts )
    else:
        id_list_noexpand = []

    print("     Building id list...")
    id_list = ci.increase_id_list( id_list_base, proj_opts, 1, 1 )
    for id in id_list_noexpand:
        id_list.append(id)
    id_list = list(set(id_list))

    ci.write_skeletons_from_list( id_list, proj_opts )

if __name__ == "__main__":
    main()
