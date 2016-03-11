function annotation_map = get_annotation_map( project, annotation_list, p )
% Generate a containers.map object mapping annotation names to skeleton
% ids.
% project: string specifying project as set up in project_parser.py.
% annotation_list: cell array of strings for annotations to query.
% p: results of python_setup_for_catmaid.m with python interpreter and
% script directories.

annotation_base = strrep(annotation_list,' ','\ ');
annotation_base = strrep(annotation_base,'''','\''');
annotation_str = '';
for ii = 1:length(annotation_base)
    annotation_str = cat(2, annotation_str, [annotation_base{ii} ' ']);
end

% This is a necessary thing to prevent MATLAB from overriding the system paths for python
oldenvdata = getenv('DYLD_LIBRARY_PATH');
setenv('DYLD_LIBRARY_PATH','')

% Run the annotation_map.py function in python with the options
eval(['!' p.py_dir ' ' p.script_dir 'annotation_map.py ' project ' ' annotation_str]);
setenv('DYLD_LIBRARY_PATH',oldenvdata);

% Read the resulting file that comes out of annotation_map.py
fid = fopen('./annotation_ids.txt');

% Walk throught the annotation_ids file and turn it into a containers.Map file for easy look-up of
% skeleton ids associated with a given annotation string.
while ~feof(fid)
    ltxt = fgetl(fid);
    cloc = regexp(ltxt,',','once');
    anno_txt = ltxt(1:cloc-1);
    anno_id_str = strsplit(ltxt(cloc+1:end),',');
    anno_id = zeros(size(anno_id_str));
    for ii = 1:length(anno_id_str);
        anno_id(ii) = str2num(anno_id_str{ii});
    end
    if ~exist('annotation_map')
        annotation_map = containers.Map({anno_txt},{anno_id});
    else
        annotation_map = [annotation_map; containers.Map({anno_txt},{anno_id})];
    end
end
