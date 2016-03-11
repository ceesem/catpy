function all_data = import_catmaid_by_annotation_v2(annotation_base,annotation_to_expand,tag_list,project,p,no_dl_flag)
% Master function for importing data from CATMAID into matlab via a folder full of JSON
% files.
% annotation_base: a cell array of strings containing annotations. Neurons
% with these annotations will be imported.
% annotations_to_expand: a cell array of strings containing annotations.
% Neurons with these annotations as well as their synaptic partners will be
% imported.
% tag_list: a cell array of strings containing node tags. The locations of
% these tags will be included in the imported data.
% project: the output of python_setup_for_catmaid.m, a struct containing the
% locations of the python interpreter and the python scripts used here.
% This file will need to be customized before this importer can run.
% no_dl_flag: Boolean value, defaults to 0 without user definition. Set to
% 1 to avoid downloading and import a folder full of json files.


% Check to see if the no_dl flag is there, set it to false if not.
if nargin == 5
    no_dl_flag = 0;
elseif no_dl_flag ~= 1
    no_dl_flag = 0;
end

if no_dl_flag == 0
    % Make a new folder for the files we're about to download.
    temp_folder = 'temp_for_catmaid';
    curr_folder = pwd();
    [~,response] = mkdir(temp_folder);
    if ~isempty(response)
        eval(['!rm ' temp_folder '/*.*'])
    end
    cd(temp_folder);
    
    %Escape out any spaces or apostrophes in the annotations
    annotation_base_es = strrep(annotation_base,' ','\ ');
    annotation_to_expand_es = strrep(annotation_to_expand,' ','\ ');
    annotation_base_es = strrep(annotation_base_es,'''','\''');
    annotation_to_expand_es = strrep(annotation_to_expand_es,'''','\''');

    
    %The script we are about to takes argument pairs that start with an
    %annotation string and then have a 0/1 flag that tells false/true for
    %expanding the download list to synaptic partners.
    
    annotation_str = '';
    for ii = 1:length(annotation_base_es)
        annotation_str = cat(2, annotation_str, [annotation_base_es{ii} ' 0 ']);
    end
    
    for ii = 1:length(annotation_to_expand)
        annotation_str = cat(2, annotation_str, [annotation_to_expand_es{ii} ' 1 ']);
    end
    
    %Run the python script after dealing with a stupid matlab library problem.
    oldenvdata = getenv('DYLD_LIBRARY_PATH');
    setenv('DYLD_LIBRARY_PATH','')
    eval(['!' p.py_dir ' ' p.script_dir 'skeletons_from_annotation.py ' project ' ' annotation_str]);
    setenv('DYLD_LIBRARY_PATH',oldenvdata);
end

% Make a list of the directory's files and initialize a dict-like Map for connector ids.
files = dir();
conn_map = containers.Map(-1, {[],[]});

disp('Importing skeletons...')

% Figure out what files are actually skeletons, since they will start with
% a given bit of text.
sk_inds = find(strncmp('sk_',{files.name},3));

% Import each json file into a struct
for ii = 1:length(sk_inds)
    disp(['    ' num2str(ii) ' of ' num2str(length(sk_inds)) '...'])
    
    % Load the json file into a cell array (of cell arrays...)
    skeleton_data = fromjson( fileread( files(sk_inds(ii)).name ) );
    
    % Turn the json file into an honest neuron, while also keeping up the
    % big list of connectors that we'll need later.
    
    [neuron, conn_map] = import_json_neuron_fast(skeleton_data, tag_list, conn_map, ii);

    
    if ~exist('all_data','var')
        all_data = neuron;
        all_data(length(sk_inds)) = neuron;           %Preallocate structure type
    else
        all_data(ii) = neuron;
    end
end

disp('Populating synapses...')
% Each neuron has to know what other neurons it talks to, so we need to go
% back through the list once more to connect them together, now that we
% have all the connectors in a containers.Map.
for ii = 1:length(all_data)
    if ~isempty(all_data(ii).synsout.connind)
        pre_conns = values(conn_map,num2cell(all_data(ii).synsout.connind));
        for jj = 1:length(pre_conns)
            all_data(ii).synsout.targinds{jj} = pre_conns{jj}{2};
        end
    end
    
    if ~isempty( all_data(ii).synsin.connind )
        
        post_conns = values(conn_map,num2cell(all_data(ii).synsin.connind));
        for jj = 1:length(post_conns)
            all_data(ii).synsin.origind(jj) = post_conns{jj}{1};
        end
        
    end
end

% This part will figure out what neurons belong to what annotations.

disp('Cleaning up skeletons...')

if no_dl_flag == 0
    cd(curr_folder)
    eval(['!rm -r ' temp_folder])
end


end

