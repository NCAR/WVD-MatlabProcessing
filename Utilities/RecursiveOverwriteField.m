% Written by: Robert Stillwell 
% Written for: NCAR
% Modificication Info: Created May 17, 2021

function [Struct] = RecursiveOverwriteField(Struct,FN,Data,Offset)
%
% Inputs: Struct: A structure to search through for data to be overwritten
%         FN:     A field name to look for for potential overwritting
%         Data:   Data to overwrite the current field with
%         Offset: An offset value to add (used for wavemeter in this case)
%
% Outputs: Struct: A structure that has been searched with possible data
%                  overwritting
%
%% Overwritting data
% Determining the fieldnames of the structure
F = fieldnames(Struct); 
% Looping over the fields to see if they are deeper structures
for m=1:1:size(F,1)
    % Checking if the current structure element is a structure
    if isstruct(Struct.(F{m,1}))
        % If a structure...dive deeper
        Struct.(F{m}) = RecursiveOverwriteField(Struct.(F{m}),FN,Data,Offset);
    else
        % Not a structure so check if 1) Is the current field the correct
        % name and 2) Is the data size the same for replacing or if that is
        % not true, check if a constant offset is to be applied
        if not(isempty(Data))
            if strcmp(F{m,1},FN) && all(size(Struct.(F{m,1})) == size(Data))
                Struct.(F{m,1}) = Data;
            end
        elseif not(isempty(Offset))
            if strcmp(F{m,1},FN)
                Struct.(F{m,1}) = Struct.(F{m,1}) + Offset;
            end
        end
    end
end
end