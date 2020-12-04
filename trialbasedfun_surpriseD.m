function [cfgin] = trialbasedfun_surpriseD(cfgin)
% TRIALBASEDFUN_SURPRISED trial definition function
%   cfgin.trialfun = 'trialfun_surpriseD'
%   for trial based analysis

% Path behavioral data
dir_behav = '/mnt/homes/home024/btalluri/confirmation_spatial/data/behavior/';

% Pull out important information about dataset to load
session = cfgin.session;
fileNum = cfgin.fileNum;
ID = cfgin.ID;

if strcmp(cfgin.subj, 'Pilot03')
    subject = 'dataRvB';
elseif strcmp(cfgin.subj, 'Pilot04')
    subject = 'dataJCT';
elseif strcmp(cfgin.subj, 'Pilot06')
    subject = 'dataAA';
end
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
idx_file = find(strcmp(files, cfgin.dataset));


[~,drug_H] = xlsread('/mnt/homes/home024/jschipp/Surprise_Drug/drug_labelAdj.xlsx');
subji = find(strcmp(ID(1:3),drug_H(:,1)));
drug_cond = drug_H(subji,2);

%%
% =====================================================================
%   INFORMATION FROM BEHAVIORAL DATA
% =====================================================================

fprintf('\n\n ---------------- \n Processing behavioral information...\n ---------------- \n\n');

% Check the run
firstblock = info_EL_blocks(idx_file,1);
lastblock = info_EL_blocks(idx_file,2);
nblock = info_EL_blocks(idx_file,1)-1;

% numFile = help variable to deal with exceptions; "fileNum" must remain
% unchanged for loading the correct file
numFile = fileNum;

% Define specific path
dir_behaviour = [dir_behav subject '/S' session+1 '/Behaviour/'];
dir_samples = [dir_behav subject '/S' session+1 '/Sample_seqs/'];

firstblock = firstblock + 1;
lastblock = lastblock-2;

trials_info = [];

% trial info to save: trial start, trial end, int1 start-500 ms, int1 end+500 ms, int2
% start-500 ms, int2 end+500 ms, binary choice or attention cue, sample mean int2, gen
% mean trial, binary choice, binary choice accuracy, estimation response.
% Petzer behav: [stim resp ACC RT start_time sacc];
% BCT behav: {'trial_start', 'trial_start_int1', 'trial_start_int2', 'gen_mean_trial', 'sample_mean_trial', 'sample_mean_int1', 'sample_mean_int2', 'interm_resp', 'interm_resp_accu', 'interm_rt', 'interm_cue', 'estim_resp', 'estim_resp_accu', 'estim_rt', 'num_sacc'};

for block = firstblock:lastblock % Loop through files/blocks
    if session < 3
        load([dir_behaviour subject(5:end) '_' session+1 '_estimation_attention_' num2str(block) '.mat'])
        load([dir_samples subject(5:end) '_' session+1 '_estimation_attention_' num2str(block) '.mat'])
    else
        load([dir_behaviour subject(5:end) '_' session+1 '_estimation_intermResp_' num2str(block) '.mat'])
        load([dir_samples subject(5:end) '_' session+1 '_estimation_intermResp_' num2str(block) '.mat'])
    end
    for trial = 1:length(Behav) % Loop through trials | Pull out info
        numSample = 12; % fixed number of samples in every block
        new_trl = [block trial numSample 0 Behav(trial, 1:4)]; % Store info
        trials_info = [trials_info; new_trl]; % Conactenate trials
    end
end



%%
% =====================================================================
%   INFORMATION FROM MEG DATA
% =====================================================================

% _____________________________________________________________________
%   TRIGGERS OF INTEREST

trg.trialOn = 11;               % Onset of pre-sequence mask
trg.respCue = 41;               % Cue for response
trg.resp = [42 43 44];          % Left/right/bad response
trg.fb = [51 52 53];            % Correct/error/bad response
trg.end = 61;                   % Offset of feedback/onset of break period
trg.sampleOn = 21;              % Onset individual sample
trg.startBlock = 1;             % Start of block
trg.endBlock = 2;               % End of block

% Variable for trigger indices
trg_idx.trialOn = [];
trg_idx.respCue = [];
trg_idx.resp = [];
trg_idx.fb = [];
trg_idx.end = [];
trg_idx.sampleOn = [];
trg_idx.startBlock = [];
trg_idx.endBlock = [];

trl = []; % function output, matrix to store trial information
% _____________________________________________________________________


% Path data - NOTE: subject tested by Julia located in a
% different folder
if cfgin.idx < 229
    cd '/mnt/homes/home024/pmurphy/meg_data/surpriseD';
else
    cd '/mnt/homes/home024/jschipp/Surprise_Drug/meg_data/';
end
% Read header and event information

hdr = ft_read_header(cfgin.dataset);
event = ft_read_event(cfgin.dataset);

% Search for trigger events - UPPT001 = trigger channel
trgVal = [event(find(strcmp('UPPT001',{event.type}))).value]';
trgIdx = find(strcmp('UPPT001',{event.type}))';

% Find events with triggers of interest and save into trg_idx structure
fields = fieldnames(trg);

for trg_type = 1:numel(fields)
    for trg_event = 1:length(trgVal)
        if ismember(trgVal(trg_event),trg.(fields{trg_type}))
            trg_idx.(fields{trg_type}) = [trg_idx.(fields{trg_type}); trgIdx(trg_event)];
        end
    end
end

% _________________________________________________________________________
% !!!   CHECK TRIGGERS & DEAL WITH EXCEPTIONS   !!!
fprintf('\n\n ---------------- \n Check triggers...\n ---------------- \n\n');

% Check whether the actual end is interrupted
if trg_idx.trialOn(end) > trg_idx.respCue(end)
    trg_idx.trialOn(end) = []; % Remove last element
end

if length(trg_idx.startBlock) > length(trg_idx.endBlock)
    trg_idx.endBlock = [trg_idx.endBlock; trg_idx.fb(end)]; % Consider the last feedback as the end of the block
end

% Check whether participant has started late
if trg_idx.trialOn(1) > trg_idx.respCue(1)
    trg_idx.respCue(1) = [];
    %trg_idx.resp(1) = [];
end

if length(trg_idx.startBlock) < length(trg_idx.endBlock)
    trg_idx.startBlock = [1; trg_idx.startBlock]; % Add a starting point
end

% Check duration of block and discard if < 3 minutes
% In some cases there are 6 blocks identified because the task has started
% without eyelink, actually only 4 of them are usable (e.g. EJG-1)
startBlockSamp = [event(trg_idx.startBlock).sample];
endBlockSamp = [event(trg_idx.endBlock).sample];
for i = 1:length(trg_idx.startBlock)
    if (endBlockSamp(i)-startBlockSamp(i))/hdr.Fs < 180 % 180 seconds
        trg_idx.startBlock(i) = 0; % Set the start/end trigger samples that are not useable to zero to remove them afterwards
        trg_idx.endBlock(i) = 0;
    end
end

% Remove start and end triggers for those blocks
startBlockSamp(trg_idx.startBlock == 0) = [];
endBlockSamp(trg_idx.endBlock == 0) = [];
trg_idx.startBlock(trg_idx.startBlock == 0) = [];
trg_idx.endBlock(trg_idx.endBlock == 0) = [];

% Remove response go cues and trial onset triggers ebfore actual block start
trg_idx.trialOn = trg_idx.trialOn(trg_idx.trialOn >= trg_idx.startBlock(1));
trg_idx.respCue = trg_idx.respCue(trg_idx.respCue >= trg_idx.startBlock(1));
%trg_idx.resp =  trg_idx.resp(trg_idx.resp >= trg_idx.startBlock(1));

% Check whether there are two response go cues between two trial onsets
% e.g. EJG-1
for i = 1:length(trg_idx.trialOn)
    if trg_idx.respCue(i)<trg_idx.trialOn(i)
        i
        trg_idx.respCue(i) = [];
        %trg_idx.resp(i) = [];
    end
end

% % Check if respCue corresponds to trialOn otherwise remove trial
% if length(trg_idx.respCue) == length(trg_idx.trialOn)
%     for i =1:length(trg_idx.respCue)
%         if (trg_idx.respCue(i) -  trg_idx.trialOn(i))>30
%             trg_idx.respCue(i) = 0;
%             trg_idx.trialOn(i) = 0;
%             i
%         end
%     end
% end
%
% trg_idx.respCue(trg_idx.respCue == 0) = [];
% trg_idx.trialOn(trg_idx.trialOn == 0) = [];



% if length(trg_idx.resp) > length(trg_idx.trialOn)
% for i = 1:length(trg_idx.trialOn)
%     if trg_idx.resp(i)<trg_idx.trialOn(i)
%         trg_idx.resp(i) = [];
%         i
%     end
% end
% end



%__________________________________________________________________________

% Initialize counters
ntrial = 0;         % Counter for trials within dataset
currBlock = 0;      % End sample of current block

fprintf('\n\n ---------------- \n Loop through blocks...\n ---------------- \n\n');
% Loop through the blocks

for start_block = trg_idx.startBlock'
    err = 0;
    block_trl = [];
    ntrial_meg = 0;
    
    nblock = nblock+1;
    
    if strcmp(ID,'UDK-1_03') && nblock == 4
        nblock = 5;
    elseif strcmp(ID,'UOC-1_01')  && nblock == 1
        nblock = nblock+2;
    elseif strcmp(ID,'UOC-1_03') && nblock == 3
        nblock = nblock+4;
    end
    
    fprintf('\n\n ---------------- \n Loop through blocks #%d\n ---------------- \n\n', nblock);
    currBlock = trg_idx.endBlock(find(trg_idx.endBlock > start_block, 1, 'first'));
    
    % Calculate the real number of trial within block based on behavioral data, starting from the end of the "trialOn" list
    % Assuming that MEG data always contains the end of recording, but sometimes the beggining is missing)
    
    % Extract blocks from behavioral data corresponding to the current file
    trials_behavBlock = trials_info(find(trials_info(:,1) == nblock),:);   % Selecting trials for this block
    
    % Create vectors containing only samples from current meg block
    trialOn_block = trg_idx.trialOn(find(trg_idx.trialOn <= currBlock & trg_idx.trialOn >= start_block));
    respCue_block = trg_idx.respCue(find(trg_idx.respCue <= currBlock & trg_idx.respCue >= start_block));
    resp_block =  trg_idx.resp(find(trg_idx.resp <= currBlock & trg_idx.resp >= start_block));
    fb_block = trg_idx.fb(find(trg_idx.fb <= currBlock & trg_idx.fb >= start_block));
    sampleOn_block = trg_idx.sampleOn(find(trg_idx.sampleOn <= currBlock & trg_idx.sampleOn >= start_block));
    
    % Remove one trialOn_block, if block was interrupted, SInce no response
    % was given in the last trial for pressing "Esc"
    if length(respCue_block) > length(resp_block)
        trialOn_block(end) = [];
        respCue_block(end) = [];
    end
    
    % !!! Exception - Something weird happend in trial 62, to make sure
    % that alignment is still correct, remove trials from 62-end
    if strcmp(ID,'XUE-2_04') && nblock == 7
        respCue_block(62:end) = []; trialOn_block(62:end) = [];
        trials_behavBlock = trials_behavBlock(1:61,:);
    end
    
    % !!! Exception - recording started after 7th trials
    if strcmp(ID,'HFK-1_02') && nblock == 1
        trials_behavBlock = trials_behavBlock(7:end,:);
    end
    
    % !!! Exception
    if strcmp(ID,'CYK-3_05') && nblock == 6
        trials_behavBlock = trials_behavBlock(11:end,:);
    end
    
    % Setting the sample from MEG data into "trial_info"
    trl_behavIdx = length(trials_behavBlock);
    
    %     if strcmp(ID,'XUE-2_04') && nblock == 7
    %         loop_vec = length(trialOn_block):-1:1;
    %         loop_vec(62) = [];
    %     else
    loop_vec = length(trialOn_block):-1:1;
    %end
    
    for i = loop_vec
        %for i =  length(trialOn_block):-1:63
        % Check whether numer of samples are the same for meg & behav
        num_samp_trl = find(sampleOn_block>trialOn_block(i) & sampleOn_block<respCue_block(i));
        if length(num_samp_trl) <= 10 % Something went wrong - more than 10 samples between trialOn and respCue
            while ~(length(num_samp_trl) ==  trials_behavBlock(trl_behavIdx,3))
                trl_behavIdx = trl_behavIdx - 1;
                if trl_behavIdx < 0, err = 1; break, end
            end
            if err == 1, break, end
            trials_behavBlock(trl_behavIdx,4) = trialOn_block(i);   % setting the sample from meg data in behav trials
        end
        trl_behavIdx = trl_behavIdx - 1;
        if trl_behavIdx < 0, err = 1; break, end
    end
    
    fprintf('\n\n ---------------- \n Loop through trials...\n ---------------- \n\n');
    % Create trials based on the trialOn trigger
    for i = trialOn_block'
        ntrial_meg = ntrial_meg + 1;
        fprintf('\n... Loop through trial %d', ntrial_meg);
        
        % Obtain trial and block info from behavioral data
        trl_behavIdx = find(trials_behavBlock(:,4) == i,4);
        %ntrial_block = trials_behavBlock(trl_behavIdx, 2);
        
        % Determine where the trial starts with respect to the event
        if isfield(cfgin.trialdef, 'prestim')
            trlOff = round(-cfgin.trialdef.prestim*hdr.Fs);
            trlStart = event(i).sample + trlOff; % Shift trial start
        else
            trlOff = event(i).offset;
            trlStart = event(i).sample;
        end
        
        % Determine where the trial ends (default response cue)
        k = respCue_block(ntrial_meg);
        if isfield(cfgin.trialdef, 'poststim')
            trlEnd = event(k).sample + round(cfgin.trialdef.poststim*hdr.Fs);  % Shift trial end
        else
            trlEnd = event(k).sample;
        end
        
        % Concatenate information
        new_trial = [trlStart trlEnd trlOff str2num(session) trials_behavBlock(ntrial_meg,:)];
        block_trl = [block_trl; new_trial];
        
    end % End of meg trial loop
    
    trl = [trl; block_trl];
    
end % End of meg block loop

trl(:,8) = []; % delete event indices

% Construct matrix containing start and end points of each block (+-10 seconds)
offset_mat = -round(10*hdr.Fs);
% Shift trial start
startSamples_mat = max([ones(length(startBlockSamp),1),(startBlockSamp + offset_mat)'],[],2);
endSamples_mat = endBlockSamp + round(10*hdr.Fs); % Shift trial end
endSamples_mat = min([event(end).sample*ones(length(endBlockSamp),1), (endBlockSamp + round(10*hdr.Fs))'],[],2);

% Extract drug condition for this session
drug_cond = drug_cond{1};
drug = drug_cond(str2double(cfgin.session));

cfgin.trl = trl;
cfgin.trialInfoLabel = {'startSample','endSample','offset','session','block', 'trial',...
    'samples','correct_choice','subjects_choice','correct_false','RT','hazard_rate','drug'};
cfgin.blockBound_trl = [startSamples_mat endSamples_mat offset_mat*ones(length(startSamples_mat),1)];
cfgin.drug = drug;
cfgin.HR = gen.H;
cfgin.event = event;

end % End of entire function

