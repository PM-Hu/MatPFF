function filedir = mkResultsDir(modelName)
% modelName - Create file folder [modelname]
% % ** code by P.M.H @bit.edu.cn (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com

% file path operation
filepath = mfilename('fullpath'); % file path
floc = strfind(filepath, '\'); % find \ 
mainpath = filepath(1 :floc(end));
filedir = fullfile([mainpath, modelName, '\']);
mkdir(filedir);

end