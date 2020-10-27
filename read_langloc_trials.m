function trl = read_langloc_trials(cfg)

%%% For langloc data %%%

%%% This function reads ECog information from .dat files      %%%
%%% INPUT: cfg = configuration file as used for fieldTrip     %%%
%%% OUTPU: trl = Nx3 matrix, where N is the number of trials  %%%
%%%              col1 = beginning of trial in the dataset     %%%
%%%              col2 = end of trial in the dataset           %%%
%%%              col3 = trigger offest within the trial       %%%
%%%               (0 = first sample of trial is the trigger)  %%%
%%%               (>0 = first sample is later than trigger)   %%%
%%%               (<0 = first sample is earlier than trigger) %%%

[ ~, states, parameters ] = load_bcidat(cfg.dataset);
sRate = parameters.SamplingRate.NumericValue;   % sampling rate
stimuli_squence=parameters.Sequence.NumericValue;
trials_value=parameters.Stimuli.NumericValue;
stimuli_value=parameters.Stimuli.Value;
% step 1: find start and end of trials
trials_indx=cell2mat(cellfun(@(x) strcmp(x,'TrialNumber'),parameters.Stimuli.RowLabels,'UniformOutput',false));
caption_indx=cell2mat(cellfun(@(x) strcmp(x,'caption'),parameters.Stimuli.RowLabels,'UniformOutput',false));
wordtype_indx=cell2mat(cellfun(@(x) strcmp(x,'Condition'),parameters.Stimuli.RowLabels,'UniformOutput',false));
StimType_indx=cell2mat(cellfun(@(x) strcmp(x,'StimType'),parameters.Stimuli.RowLabels,'UniformOutput',false));
IsRight_indx=cell2mat(cellfun(@(x) strcmp(x,'IsProbeCorrect'),parameters.Stimuli.RowLabels,'UniformOutput',false));
trial_for_stimuli_seq=trials_value(trials_indx,:);
trials=unique(trial_for_stimuli_seq);
trials(isnan(trials))=[];
fprintf('%d trials were found \n',length(trials));
trial_seq_cell={};
for ii=1:length(trials)
    trial_stimuli_sequence=find(trial_for_stimuli_seq==trials(ii));
    trial_instance_in_sequence=strfind(stimuli_squence',trial_stimuli_sequence);
    if length(trial_instance_in_sequence)==1
        trial_seq_cell{ii,1}=stimuli_squence(trial_instance_in_sequence+[0:length(trial_stimuli_sequence)-1]);
        trial_seq_cell{ii,2}=trials(ii);
    elseif isempty(trial_instance_in_sequence)
        fprintf('the stimuli for trial %d was not find in the parameter.sequence \n',trials(ii))
    else
        fprintf('more than on instance of trial found\n');
        keyboard;
    end
end
trl = zeros(length(length(trial_seq_cell)),3);
for k=1:length(trial_seq_cell)
    trial_indx=trial_seq_cell{k};
    % find trial type
    wordtype=trials_value(find(wordtype_indx),trial_for_stimuli_seq==trial_seq_cell{k,2});
    stimuli_range=[];
    fprintf('adding trial %d  \n',(k))
    for kk=1:length(trial_indx)
        stimulus_index=find(states.StimulusCode==trial_indx(kk));
        stimuli_range=[stimuli_range;[min(stimulus_index),max(stimulus_index)]];
    end
    trl(k,1)=min(stimuli_range(:))-sRate*cfg.trialdef.prestim;
    trl(k,2)=max(stimuli_range(:))+sRate*cfg.trialdef.poststim;
    trl(k,3)=-sRate*cfg.trialdef.prestim-1;
end



end