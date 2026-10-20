function oldStruct = UpdateStruct(oldStruct, newStruct)
%% Update a structure using a new structure provided. This new structure could be a subset of the old structure but must have the same fieldnames

    oldNames = fieldnames(oldStruct) ; 
    newNames = fieldnames(newStruct) ; 
    for n=1:length(newNames)
        % ensure structure fieldnames match
        if ~ismember(newNames{n}, oldNames) 
            error('%s: The new structure contains a field that is not in the old structure', mfilename) ; 
        end
        % 
        if ~isstruct(newStruct.(newNames{n}))
            oldStruct.(newNames{n}) = newStruct.(newNames{n}) ; 
        else % the field is a structure so we need to go one level deeper
            oldStruct.(newNames{n}) = mati.Tools.UpdateStruct(oldStruct.(newNames{n}), newStruct.(newNames{n})) ; 
        end
        
    end




%% 
