function [neuron, conn_map] = import_json_neuron_fast( data, tag_list, conn_map, sk_ind )
% Import a neuron into matlab from its JSON as extracted with the
% 'fromjson.m' function in matlab-json package.

poststruct = struct('xyz',zeros(0,3),'orig_id',[],'origind',[],'treeinds',[],'connind',[]);
prestruct = struct('xyz',zeros(0,3),'numTargs',[],'treeinds',[],'connind',[]);

% Massage the node data into a proper matrix
node_cat = cat(2, data{1}{:});
rootnodeind = find(cellfun(@(x)isnan(x),node_cat(2,:)));
node_cat{2,rootnodeind(1)} = -1;
nodemat = [ [node_cat{1,:}];
    [node_cat{2,:}];
    [node_cat{3,:}];
    [node_cat{4,:}];
    [node_cat{5,:}];
    [node_cat{6,:}];
    [node_cat{7,:}];
    [node_cat{8,:}];
    ]';

neuron = struct('A',sparse(size(nodemat,1)),'Aw',sparse(size(nodemat,1)),'xyz',double( nodemat(:,4:6) ),'R',nodemat(:,7),'name',data{5},'L',[],'root',rootnodeind,'id',data{4},'soma',-1,'vids',nodemat(:,1),'synsin',poststruct,'synsout',prestruct);
neuron.synsout.targinds = cell(0,1);
neuron.synsout.targ_ids = cell(0,1);

nonrootinds = [1:rootnodeind-1 rootnodeind+1:size(nodemat,1)];

map = containers.Map(nodemat(:,1),1:length(nodemat(:,1)));
parents = arrayfun(@(x)map(x),nodemat(nonrootinds(:),2));

A = sparse( nonrootinds, parents', ones(size(parents')), size(nodemat,1), size(nodemat,1));

neuron.A = A + A';
neuron.Aw = spatializeAdjMat(neuron);
neuron.L = sum(neuron.Aw(:))/2;
neuron.Adir = sparse( size(A,1), size(A,1));
neuron.Adir( A==1 ) = neuron.Aw( A==1 );
neuron.id2ind = map;

%     # 1: treenode_connector.treenode_id
%     # 2: treenode_connector.connector_id
%     # 3: 0 for presynaptic, 1 for postsynaptic
%     # 4: connector.location.x
%     # 5: connector.location.y
%     # 6: connector.location.z

if ~isempty(data{2})
    conn_grid = cat(2,data{2}{:});
    conn_mat = [ [conn_grid{1,:}];...
        [conn_grid{2,:}];...
        [conn_grid{3,:}];...
        [conn_grid{4,:}];...
        [conn_grid{5,:}];...
        [conn_grid{6,:}]...
        ]';
    
    
    connperneuron = conn_mat(:,2);
    isoldkey = isKey(conn_map,num2cell(connperneuron));
    
    post_inds = conn_mat(:,3) == 1;
    neuron.synsin.connind = conn_mat(post_inds,2);
    neuron.synsin.treeinds = arrayfun(@(x)map(x),conn_mat(post_inds,1));
    neuron.synsin.xyz = double( conn_mat(post_inds,4:6) );
    
    for ii = 1:length(neuron.synsin.connind)
        if isKey(conn_map,neuron.synsin.connind(ii))
           base_conn = conn_map(neuron.synsin.connind(ii));
        else
           base_conn = {-1,[]};
        end
        base_conn{2}(end+1) = sk_ind;
        conn_map(neuron.synsin.connind(ii)) = base_conn;
        %conn_struct(conn_map(neuron.synsin.connind(ii))).post(end+1) = sk_ind;
    end
    
    pre_inds = conn_mat(:,3) == 0;
    neuron.synsout.connind = conn_mat(pre_inds,2);
    
    neuron.synsout.treeinds = arrayfun(@(x)map(x),conn_mat(pre_inds,1));
    neuron.synsout.xyz = double( conn_mat(pre_inds,4:6) );
    
    if ~isempty(data{6})
        synweight_grid = cat(2,data{6}{:});
        synweightmat = [[synweight_grid{1,:}];...
            [synweight_grid{2,:}]]';
        local_conn_map = containers.Map(neuron.synsout.connind,1:length(neuron.synsout.connind));
        inds_for_weight = arrayfun(@(x)local_conn_map(x),synweightmat(:,1));
        neuron.synsout.numTargs(inds_for_weight) = synweightmat(:,2)';
    end
    
    for ii = 1:length(neuron.synsout.connind)
        if isKey(conn_map,neuron.synsout.connind(ii))
           base_conn = conn_map(neuron.synsout.connind(ii));
           base_conn{1} = sk_ind;
        else
           base_conn = {sk_ind,[]};
        end
        conn_map(neuron.synsout.connind(ii)) = base_conn;

%         conn_struct(conn_map(neuron.synsout.connind(ii))).pre = sk_ind;
    end
end

% Note the location of specific tags that have been given to the neuron.
% For a few of the hardwired tags, add them by default.
tag_list = unique([tag_list {'soma','ends','not a branch','microtubules end'}]);
neuron = return_tag_locations( neuron, data, tag_list);
