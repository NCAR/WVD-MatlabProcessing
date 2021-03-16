% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created August 18th, 2020

function [IsField,Data] = RecursivelyCheckIsField(Struct, Field)
%
% Input: Struct:    A structure to check for the existance of a field
%        Field:     Either a string or a cell array of strings to check if
%                   a field exists in a given structure
%
% Outputs: IsField: Boolean value specifying if the requested field exists
%          Data:    Copy of the data of the field requested if availible
%
%% Assume that the field does not exist
IsField = false;
Data    = [];
%% Recursively check for the existance of the field
if isa(Field,'char') || length(Field) == 1
    % At the bottom of the cell array so check on data
    IsField = isfield(Struct,Field);
    % If data exists, make a copy to send back 
    if IsField
        try
            Data = Struct.(Field{:});
        catch
            Data = Struct.(Field);
        end
    end
else
    % Not at the end of the cell array so dive down further if possible
    if isfield(Struct,Field{1})
        [IsField,Data] = RecursivelyCheckIsField(Struct.(Field{1}),Field(2:end));
    end
end
end