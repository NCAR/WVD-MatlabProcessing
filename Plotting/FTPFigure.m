% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created July 2nd, 2021

function FTPFigure(FN,Op,Paths,Serv,Type)
%
% Inputs: FN:    The figure number to save and close
%         Op:    A structure containing all of the user defined options
%         Paths: A structure containing all the user defined file path info
%         Serv:  A boolean indicating if code is running on the server
%         Type:  A string containing the type of figure to be saved
%
%%
if Op.UploadFig && Serv
    % Creating the name of the file to submit
    Name = ['lidar.NCAR_',upper(erase(Op.System,'_')),'.',Op.Date,'00.',Type,'.png'];
    % Printing a temp file for uploading
    print(FN, Name, '-dpng', '-r300') % set the resolution as 300 dpi
    % Opening up an FTP connection to the field catalog and uploading file
    test=sftp('catalog.eol.ucar.edu', "rsfdata", "Password", "MakeHcrLookGood5#");
    cd(test,'operations');
    mput(test, Name);
    close(test);
    % Removing the temp file so as not to clog up the working directory
    delete(Name)
end
end
