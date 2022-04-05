% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created July 2nd, 2021

function FTPFigure(FN,Op,Paths,Type)
%
% Inputs: FN:    The figure number to save and close
%         Op:    A structure containing all of the user defined options
%         Paths: A structure containing all the user defined file path info
%         Type:  A string containing the type of figure to be saved
%
%%
if Op.UploadFig
    % Creating the name of the file to submit
    Name = ['lidar.NCAR_',upper(erase(Op.System,'_')),'.',Op.Date,'00.',Type,'.png'];
    % Printing a temp file for uploading
    print(FN, Name, '-dpng', '-r300') % set the resolution as 300 dpi
    % Opening up an FTP connection to the field catalog and uploading file
    test=ftp('catalog.eol.ucar.edu', 'anonymous', 'spuler@ucar.edu');
    cd(test,Paths.FieldCat);
    mput(test, Name);
    cd(test);
    dir(test,'lidar*')
    close(test);
    % Removing the temp file so as not to clog up the working directory
    delete(Name)
end
end
