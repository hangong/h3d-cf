function list = getAll(dirName,tp)
% get a list of files in dirName
% tp - whether to obtain files

dirData = dir(dirName);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
if strcmp(tp,'f')
    idx = ~dirIndex;
else
    idx = dirIndex & ~ismember({dirData.name},{'.','..'});
end

if isempty(idx)
    list = {};
else
    list = {dirData(idx).name}';  % Get a list of the files
end
