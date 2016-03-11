# Annotation_map takes a list of strings (annotations) in the command line and
# writes a file where each row is for a different annotation.
# The first element is the annotation itself, then a comma-separated list of skeleton ids
# for the neurons having that annotation.

import catmaid_interface as ci
import cardonalab_project_parser as cpp
import sys

def main():

    project_flag = sys.argv[1]
    input_args = [ sys.argv[i] for i in range(2,len(sys.argv)) ]

    proj_opts = cpp.project_map( project_flag )
    anno_dict = ci.get_annotation_dict( proj_opts )
    arg_ids = [ anno_dict[arg] for arg in input_args ]

    f_anno = open('annotation_ids.txt','w')

    for arg in input_args:
        f_anno.write(arg)
        id_list_base = ci.get_ids_from_annotation( [anno_dict[arg]], proj_opts )
        for skid in id_list_base:
            f_anno.write(',' + str(skid))
        f_anno.write('\n')
    f_anno.close()

if __name__ == "__main__":
    main()
