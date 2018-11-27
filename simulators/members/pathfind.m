function pathfind(edgeSet, source, sink, path)
% PATHFIND is to find all paths from the given source to the given
% sink.
%
% Inputs:
%   - edgeSet: 2xK matrix, each 2-dim vector specifies an arc
%   - source: the starting node
%   - sink: the ending node
%   - path: must be set to []
% Outputs:
%   all paths are saved in the csv file "pathTable.csv"


if source == sink
    path = [path sink];
    fid = fopen('pathTable.csv', 'a');
    for k = 1:length(path)-1
        fprintf(fid, '%d,', path(k));
    end
    fprintf(fid, '%d\n', path(k+1));
    fclose(fid);
else
    idxList = find(edgeSet(1,:) == source);
    idxRest = 1:length(edgeSet(1,:));
    idxRest(idxList) = [];
    edgeSetNew = edgeSet(:,idxRest);

    path = [path source];
    for idx = idxList
        sourceNew = edgeSet(2,idx);
        pathfind(edgeSetNew, sourceNew, sink, path);
    end
end

end