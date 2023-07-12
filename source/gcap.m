function gcap(mrcc, input, output, param)
%% Processing input arguments
default_mrcc_source = 'MrCCPackage.cpp';    % source file for subspace clustering
default_mrcc_exe = 'MrCCPackage';           % exe file for subspace clustering
default_mrcc_alpha = 1e-10;                 % MrCC argument -- please check the mrcc readme
default_mrcc_h = 4;                         % MrCC argument -- please check the mrcc readme
default_mrcc_hard = 0;                      % MrCC argument -- 0/1 for soft/hard clustering
default_input_feature = '../input/feature.input';       % input file containing image features
default_input_location = '../input/location.input';     % input file containing image tile (x, y)
default_input_label = '../input/label.input';           % input file containing image labels
default_input_query = '../input/query.input';           % input file containing image queries
default_output_cluster = '../output/cluster.output';    % output file containing cluster information (as a result of the clustering algorithm)
default_output_label = '../output/label.output';        % output file containing the top label (with normalized similarity score) of each image
default_output_detail = '../output/label.detail';       % output file containing the detailed labeling result of each image
default_param_r = 1;        % r: for each image, the ratio of image->label edge weights w.r.t. image->cluster edge weights (default value is 1)
default_param_c = 0.5;      % c: (1-c) is the restart probability
default_param_l = 0;        % l: weights for geographic edges
default_param_outlier = true;   % outlier: include outlier images as an additional cluster
small = 1e-10;

if nargin < 1, 
    mrcc = struct('source', default_mrcc_source, 'exe', default_mrcc_exe, 'alpha', default_mrcc_alpha, 'h', default_mrcc_h, 'hard', default_mrcc_hard);
end;
if nargin < 2, 
    input = struct('feature', default_input_feature, 'label', default_input_label, 'query', default_input_query, 'location', default_input_location);
end;
if nargin < 3, 
    output = struct('cluster', default_output_cluster', 'label', default_output_label, 'detail', default_output_detail);
end;
if nargin < 4,
    param = struct('r', default_param_r', 'c', default_param_c, 'l', default_param_l, 'outlier', default_param_outlier);
end;

%% Perform subspace clustering
% First get the total number of features and total number of lines
nImage = 0;
fid = fopen(input.feature, 'r');
tline = fgetl(fid);
c = textscan(tline, '%f', 'delimiter', '\t');
incl_ground_truth = 1;  % Note: modify this line if ground truth is not included  
nFeature = length(c{1}) - incl_ground_truth;      
while ischar(tline)
    nImage = nImage + 1;
    tline = fgetl(fid);
end;
fclose(fid);
% Then compile and run the code
stmt = sprintf('g++ %s -o %s', mrcc.source, mrcc.exe);
unix(stmt);
stmt = sprintf('./%s %e %f %d %s %s %d %d', mrcc.exe, mrcc.alpha, mrcc.h, mrcc.hard, input.feature, output.cluster, nFeature, nImage); 
unix(stmt);
% fprintf('%s\n', stmt);

%% Graph construction
% Read the input of clustering results
fid = fopen(output.cluster, 'r');
c = textscan(fid, '%d %d', 'delimiter', ',', 'CollectOutput', 1);
c = double(c{1});
fclose(fid);
nCluster = max(c(:,2));
oindex = (c(:,2) == 0);
if param.outlier
   c(oindex, 2) = nCluster + 1;
   nCluster = nCluster + 1;
else
   c = c(~oindex, :); 
end;
image_cluster = [c(:,1) c(:,2) ones(size(c,1),1)];
% Read the input of location
fid = fopen(input.location, 'r');
w_ii = sparse(nImage, nImage);
if fid >= 0,
    c = textscan(fid, '%d %d %d', 'delimiter', '\t', 'CollectOutput', 1);
    c = c{1};
    loc = cell(max(c(:,1)), 1);
    for i = 1 : size(loc,1)
        loc{i} = sparse(0,0);
    end;
    for i = 1 : size(c, 1)
        loc{c(i,1)}(c(i,2), c(i,3)) = i;
    end;
    fclose(fid);    
    for i = 1 : size(loc, 1)
        if ~isempty(loc{i})
            [n1, n2] = size(loc{i});
            for j = 1 : n2-1
                index = find(loc{i}(:,j) .* loc{i}(:,j+1) ~= 0);
                for k = 1 : length(index)
                    i1 = loc{i}(index(k), j);
                    i2 = loc{i}(index(k), j+1);
                    w_ii(i1,i2) = 1;
                    w_ii(i2,i1) = 1;
                end;
            end;
            for j = 1 : n1-1
                index = find(loc{i}(j,:) .* loc{i}(j+1,:) ~= 0)';
                for k = 1 : length(index)
                    i1 = loc{i}(j, index(k));
                    i2 = loc{i}(j+1, index(k));
                    w_ii(i1,i2) = 1;
                    w_ii(i2,i1) = 1;
                end;
            end;
        end;
    end;
end;
% Read the input of labeling
fid = fopen(input.label, 'r');
c = textscan(fid, '%d %d', 'delimiter', ',', 'CollectOutput', 1);
c = double(c{1});
fclose(fid);
nLabel = max(c(:,2));
image_label = [c(:,1) c(:,2) ones(size(c,1),1)];
% Construct the sparse matrix
w_ci = getsparse(image_cluster, nImage, nCluster)';
w_li = getsparse(image_label, nImage, nLabel)';
for j = 1 : nImage
    w_ci(:,j) = w_ci(:,j) ./ ( sum(w_ci(:,j)) + small );
    w_li(:,j) = w_li(:,j) ./ ( sum(w_li(:,j)) + small ) .* param.r; 
    w_ii(:,j) = w_ii(:,j) ./ ( sum(w_ii(:,j)) + small ) .* param.l;
end;
[i1, j1, v1] = find(w_ci);
[i2, j2, v2] = find(w_li);
[i3, j3, v3] = find(w_ii);
nNode = nImage + nCluster + nLabel;
W = getsparse([ [i1+nImage, j1, v1]; [j1, i1+nImage, v1]; [i2+nImage+nCluster, j2, v2]; [j2, i2+nImage+nCluster, v2] ; [i3, j3, v3] ], nNode, nNode);
for j = 1 : nNode
    W(:,j) = W(:,j) ./ ( sum(W(:,j)) + small ); 
end;
%% Random Walk with Restart
fid = fopen(input.query', 'r');
c = textscan(fid, '%d');
fclose(fid);
c = c{1};
f1 = fopen(output.label, 'w');
f2 = fopen(output.detail, 'w');
fprintf('RWR Query...  0.0%%');
for i = 1 : length(c)
    r = gcap_rwr(W, param.c, c(i));
    r = r(end-nLabel+1 : end);
    [lbl, lbl] = max(r);
    fprintf(f1, '%d\t%d\t%8.6f\n', c(i), lbl, r(lbl)/sum(r));
    fprintf(f2, '%d', c(i));
    fprintf(f2, '\t%8.6f', r ./ sum(r));
    fprintf(f2, '\n');
    fprintf('\b\b\b\b\b\b%5.1f%%', 100 * i / length(c));
end;
fprintf('\n');
fclose(f1);
fclose(f2);


function r = gcap_rwr(W, c, qid)
% The graph query function for GCap
%% Input:
% W: the graph weight matrix (column normalized)
% c: 1-c is the restart probability
% qid: a vector containing the node-id of the query
%% Output:
% r: the steady-state probability vector
nNode = length(W);
e = sparse(nNode, 1);
e(qid) = 1 / length(qid);
max_iter = 100;
max_diff = 0.01 / nNode;
old_r = zeros(nNode, 1);
r = repmat(1/nNode, nNode, 1);
iter = 0;
while sqrt(mean((r - old_r) .^ 2)) >= max_diff && iter < max_iter
    iter = iter+1;
    old_r = r;
    r = c * W * r + (1-c) * e;
    r = r ./ sum(r);
end;

function w = getsparse(edgeList, nRow, nCol)
w = spconvert(double(edgeList));
w = [[w, sparse(size(w,1), nCol-size(w,2))]; sparse(nRow-size(w,1), nCol)];
