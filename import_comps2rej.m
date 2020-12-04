clear, close all

curr_folder = pwd;
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/misc/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/plotting/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/stats/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/SDT/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/fitting/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/Colormaps/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/fieldtrip_2020/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/CircStat2012a/
addpath /mnt/homes/home024/btalluri/Tools/MATLAB/eye/
addpath(curr_folder)

ft_defaults

% Load important information about files
% files contains the complete names of all files that must be processed
% info_EL_blocks: matrix with 5 columns
% col 1 - first block to process with regard to behavioral data (normally 1 or 5)
% col 2 - last block to process with regard to behavioral data (normally 4 or 8)
% col 3 - number of blocks in file (normally 4)
% col 4 - run (1 or 2)
% col 5 - eyelink data (1: usable, 0: use veog instead)
[~,files] = xlsread('Info_filewise');
info_EL_blocks = xlsread('Info_filewise');

[~,txt,~] = xlsread('Manual_Comp_Rej.xlsx');
txt = txt(:,2:4);

for i = [1:8 10:12]
    filein = files{i};
    
    ID = [filein(1:9) filein(end-5:end-3) '_' int2str(info_EL_blocks(i,4))]; % Subject ID + Session number + file number
    if regexp(ID, 'Pilot03_C*') % With this subject the session number was missing when registered
        ID = ['Pilot03_1' filein(end-5:end-3) '_' int2str(info_EL_blocks(i,4))];
    end
    
    rej = txt(i,:);
    rej = [rej{1,1} rej{1,2} rej{1,3}];
    rejComps = [];
    
    remain = rej;
    while ~isempty(remain)
        [token,remain] = strtok(remain, '; ');
        token = str2double(token);
        rejComps = [rejComps; token];
    end
    rejComps = unique(rejComps)';
    rejComps(isnan(rejComps)) = [];
    
    cd(['/mnt/homes/home024/btalluri/confirmation_spatial/data/meg/analysis/comp_ICA/' ID(1:7)])
    comp2rej = ['comp2rej_' ID '.mat'];
    save(comp2rej,'rejComps')
end
