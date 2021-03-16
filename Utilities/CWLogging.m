



function CWLogging(String,Options,Type)
%
%
%
%
%
%%
if strcmp(Options.Logging,'Full') || strcmp(Type,'Warning')
    fprintf(String)
elseif strcmp(Options.Logging,'Skinny')
    if strcmp(Type,'Main') 
        fprintf(['   ',String])
    end
end
end