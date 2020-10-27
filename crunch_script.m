%% extract dat files 
clear all 
close all 
home

%% 
exp_name='MITLangloc';

%% 
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end 
    
%% 
data_path='/Users/eghbalhosseini/MyData/ecog-langloc';
save_path='/Users/eghbalhosseini/MyData/ecog-langloc/crunched/';
if ~exist(save_path)
    mkdir(save_path);
end
d= dir([data_path,'/**/ECOG*.dat']);
fprintf(' %d .dat files were found \n', length(d))
%% 
d_op_info= dir([data_path,'/**/*operation_info.mat']);
fprintf(' %d op_info files were found \n', length(d));

%% 
for i=8%:5%length(d)
    fprintf('extracting %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subject_name=d(i).folder(strfind(d(i).folder,'AMC')+[0:5]);
    session_name=d(i).name(1:end-4);
    
    subject_op_id=find(~cellfun(@isempty,(cellfun(@(x) strfind(x,subject_name),{d_op_info(:).name},'UniformOutput',false))));
    subject_session_info=[];
    subject_op_info=[];
    exp_info=[];
    if ~isempty(subject_op_id)
        subject_op_info=load(strcat(d_op_info(subject_op_id).folder,'/',d_op_info(subject_op_id).name));
        % test whether its a good session 
        subject_op_info=subject_op_info.(strcat(subject_name,'_op'));
        exp_info=subject_op_info.exp_info.(exp_name);
        session_info=exp_info(~cellfun(@isempty,cellfun(@(x) strfind(x,session_name),exp_info.sessions,'UniformOutput',false)),:);
        is_good_session=session_info.good_run;
        subject_session_info.op_info=subject_op_info.op_info;
        subject_session_info.session_info=session_info;
    else 
        is_good_session=1;
    end
    if is_good_session
    % find the 
    subject_op_info=[];
    [signal_broadband,...
     signal_bandpass,...
     signal_envelope,...
     signal_evelope_downsample,...
     signal_hilbert_downsample,...
     signal_hilbert_zs_downsample,...
     states,...
     parameters,...
     ecog_parameters]=filter_channels_using_schalk({strcat(d(i).folder,'/',d(i).name)},subject_op_info);

    
    % start with an empty strcuture for data and info 
    list_var_to_get={''};
    stim_types={'Jabberwocky','Sentences'};
    dat={};
    info=struct;
    
    % step 1: find start and end of trials 
    info.sample_rate=parameters.SamplingRate.NumericValue;
    stimuli_squence=parameters.Sequence.NumericValue;
    trials_value=parameters.Stimuli.NumericValue;
    stimuli_value=parameters.Stimuli.Value;
    %
    trials_indx=cell2mat(cellfun(@(x) strcmp(x,'TrialNumber'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    caption_indx=cell2mat(cellfun(@(x) strcmp(x,'caption'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    wordtype_indx=cell2mat(cellfun(@(x) strcmp(x,'Condition'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    StimType_indx=cell2mat(cellfun(@(x) strcmp(x,'StimType'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    IsRight_indx=cell2mat(cellfun(@(x) strcmp(x,'IsProbeCorrect'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    %
    trial_for_stimuli_seq=trials_value(trials_indx,:);
    trials=unique(trial_for_stimuli_seq);
    trials(isnan(trials))=[];
    %fprintf('%d trials were found \n',length(trials));
    % 
    figure;
    plot(states.StimulusCode);
    %keyboard
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
    % extracting data per trial for subject
    trial_reponse=[];
    % do a transpose 
    signal_broadband=transpose(signal_broadband);
    signal_hilbert_downsample=transpose(signal_hilbert_downsample);
    signal_hilbert_zs_downsample=transpose(signal_hilbert_zs_downsample);
    
    signal_bandpass=cellfun(@transpose,signal_bandpass,'UniformOutput',false);
    signal_envelope=cellfun(@transpose,signal_envelope,'UniformOutput',false);
    signal_evelope_downsample=cellfun(@transpose,signal_evelope_downsample,'UniformOutput',false);
    
    for k=1:length(trial_seq_cell)
        trial_indx=trial_seq_cell{k};
        % find trial type 
        wordtype=trials_value(find(wordtype_indx),trial_for_stimuli_seq==trial_seq_cell{k,2});
        wordtype(isnan(wordtype))=[];
        if ~isempty(wordtype)
            info.word_type{k,1}=stim_types{unique(wordtype)};
        else
            info.word_type{k,1}='0';
        end 
        trial=struct;
        trial_index=[];
        trial_downsample_index=[];
        trial_string=[];
        trial_probe=[];
        trial_type=[];
        stimuli_range=[];
        stimuli_downsample_range=[];
        stimuli_type={};
        probe_result=[];
        stimuli_string={};
        signal_bandpass_parsed={};
        LeftEyeGazeX={};
        LeftEyeGazeY={};
        LeftEyePosX={};
        LeftEyePosY={};
        LeftPupilSize={};
        signal_broadband_parsed={};
        signal_bandpass_envelope_parsed={};
        signal_hilbert_downsample_parsed={};
        signal_hilbert_zs_downsample_parsed={};
        signal_bandpass_envelope_downsample_parsed={};
        
        fprintf('adding trial %d  \n',(k))
        for kk=1:length(trial_indx)
            stimulus_index=find(states.StimulusCode==trial_indx(kk));
            stimuli_downsample_index=find(states.StimulusCodeDownsample==trial_indx(kk));
            stimuli_type{kk,1}=stimuli_value{StimType_indx,trial_indx(kk)};
            if ~isempty(stimuli_value{IsRight_indx,trial_indx(kk)})
                probe_result=[probe_result,stimuli_value{IsRight_indx,trial_indx(kk)}];
            end 
            stimuli_string{kk,1}=stimuli_value{caption_indx,trial_indx(kk)};
            trial_index=[trial_index;stimulus_index];
            trial_downsample_index=[trial_downsample_index;stimuli_downsample_index];
            % 
            stimuli_range=[stimuli_range;[min(stimulus_index),max(stimulus_index)]];
            stimuli_downsample_range=[stimuli_downsample_range;[min(stimuli_downsample_index),max(stimuli_downsample_index)]];
            %
            signal_broadband_parsed{kk,1}=signal_broadband(:,stimulus_index);
            signal_hilbert_downsample_parsed{kk,1}=signal_hilbert_downsample(:,stimuli_downsample_index);
            signal_hilbert_zs_downsample_parsed{kk,1}=signal_hilbert_zs_downsample(:,stimuli_downsample_index);
            % 
            signal_bandpass_parsed=[signal_bandpass_parsed;transpose(cellfun(@(x) x(:,stimulus_index),signal_bandpass,'UniformOutput',false ))];
            signal_bandpass_envelope_parsed=[signal_bandpass_envelope_parsed;transpose(cellfun(@(x) x(:,stimulus_index),signal_envelope,'UniformOutput',false ))];
            signal_bandpass_envelope_downsample_parsed=[signal_bandpass_envelope_downsample_parsed;...
                transpose(cellfun(@(x) x(:,stimuli_downsample_index),signal_evelope_downsample,'UniformOutput',false ))];
            % 
            LeftEyeGazeX=[LeftEyeGazeX;transpose(states.EyetrackerLeftEyeGazeX(stimulus_index,:))];
            LeftEyeGazeY=[LeftEyeGazeY;transpose(states.EyetrackerLeftEyeGazeY(stimulus_index,:))];
            LeftEyePosX=[LeftEyePosX;transpose(states.EyetrackerLeftEyePosX(stimulus_index,:))];
            LeftEyePosY=[LeftEyePosY;transpose(states.EyetrackerLeftEyePosY(stimulus_index,:))];
            LeftPupilSize=[LeftPupilSize;transpose(states.EyetrackerLeftPupilSize(stimulus_index,:))];
            % 
           
            if strfind(stimuli_value{StimType_indx,trial_indx(kk)},'word')
                trial_string=[trial_string,' ',stimuli_value{caption_indx,trial_indx(kk)}];
            end 
            if strfind(stimuli_value{StimType_indx,trial_indx(kk)},'probe')
                trial_probe=[trial_probe,' ',stimuli_value{caption_indx,trial_indx(kk)}];
            end
            trial_type=[trial_type,' ',stimuli_value{StimType_indx,trial_indx(kk)}];
        end 
        %trial.(strcat('signal','_broadband'))=signal_broadband(trial_index,:);
        %trial.(strcat('signal','_bandpass'))=cellfun(@(x) x(trial_index,:),signal_bandpass,'UniformOutput',false );
        %trial.(strcat('signal','_envelope'))=cellfun(@(x) x(trial_index,:),signal_envelope,'UniformOutput',false );
        %trial.(strcat('signal','_envelope_downsample'))=cellfun(@(x) x(trial_downsample_index,:),signal_evelope_downsample,'UniformOutput',false );
        %trial.(strcat('signal','_broadband_parsed'))=signal_broadband_parsed;
        %trial.(strcat('signal','_bandpass_parsed'))=signal_bandpass_parsed;
        trial.(strcat('signal','_hilbert_downsample_parsed'))=signal_hilbert_downsample_parsed;
        trial.(strcat('signal','_hilbert_zs_downsample_parsed'))=signal_hilbert_zs_downsample_parsed;
        %trial.(strcat('signal','_envelope_parsed'))=signal_bandpass_envelope_parsed;
        trial.(strcat('signal','_envelope_dowsample_parsed'))=signal_bandpass_envelope_downsample_parsed;
        trial.(strcat('signal_ave','_envelope_downsample_parsed'))=cellfun(@(x) nanmean(x,2),signal_bandpass_envelope_downsample_parsed,'UniformOutput',false);
        trial.(strcat('signal_ave','_hilbert_downsample_parsed'))=cellfun(@(x) nanmean(x,2),signal_hilbert_downsample_parsed,'UniformOutput',false);
        trial.(strcat('signal_ave','_hilbert_zs_downsample_parsed'))=cellfun(@(x) nanmean(x,2),signal_hilbert_zs_downsample_parsed,'UniformOutput',false);

        %trial.(strcat('eye','_LeftEyeGazeX_parsed'))=LeftEyeGazeX;
        %trial.(strcat('eye','_LeftEyeGazeY_parsed'))=LeftEyeGazeY;
        %trial.(strcat('eye','_LeftEyePosX_parsed'))=LeftEyePosX;
        %trial.(strcat('eye','_LeftEyePosY_parsed'))=LeftEyePosY;
        %trial.(strcat('eye','_LeftPupilSize_parsed'))=LeftPupilSize;
        trial.(strcat('trial','_string'))=trial_string;
        trial.(strcat('trial','_probe_question'))=trial_probe;
        trial.(strcat('trial','_probe_answer'))=probe_result;
        trial.trial_onset_sec=stimulus_index(1)./info.sample_rate;
        %trial.keydown=states.KeyDown(trial_index);
        %trial.keyup=states.KeyUp(trial_index);
        %trial.isRight=states.IsRight(trial_index);
        trial.signal_range=stimuli_range-min(stimuli_range(:))+1;
        trial.signal_range_downsample=stimuli_downsample_range-min(stimuli_downsample_range(:))+1;
        trial.stimuli_type=stimuli_type;
        trial.stimuli_string=stimuli_string;
        % find subject response: 
        index_isright_start = trial_index(find(diff(double(states.IsRight(trial_index) > 0)) == 1)+1);
        index_isright_stop  = trial_index(find(diff(double(states.IsRight(trial_index) > 0)) == -1));
        buffer_before = info.sample_rate * 1; % 1 sec 
        buffer_after  = info.sample_rate * 2; % 2 sec
        KeyDown = unique(states.KeyDown((index_isright_start-buffer_before):(index_isright_stop+buffer_after)));
        KeyDown = intersect(KeyDown,[67,77,99,109]);
        if length(KeyDown) ~= 1                    % too many key's pressed or incorrect response
           TrialResponse = 'INCORRECT_KEY';
        elseif KeyDown == 67 || KeyDown == 99      % response is yes (1)
           TrialResponse = 'RIGHT'; 
        elseif KeyDown == 77 || KeyDown == 109     % response is no  (2) 
           TrialResponse = 'WRONG'; 
        else                                       % incorrect response 
           TrialResponse = 'INCORRECT_KEY';
        end
        % 
        trial.subject_response=TrialResponse;
        info.subject_response{k,1}=TrialResponse;
        info.probe_value{k,1}=probe_result;
        dat{k,1}=trial;
        if ~contains(trial_type,'word')
            info.trial_type{k,1}='fixation';
        else
            info.trial_type{k,1}='word';
        end
        
    end
    info.subject=subject_name;
    info.session_name=session_name;
    %
    info.random_stim_present=parameters.SequenceType.NumericValue;
    info.random_stim_present_comment=parameters.SequenceType.Comment;
    %
    info.num_of_stim_rep=parameters.NumberOfSequences.NumericValue;
    info.num_of_stim_rep_comment=parameters.NumberOfSequences.Comment;
    %
    try
        info.common_refs=parameters.CommonReference.NumericValue;
        info.common_refs_comment=parameters.CommonReference.Comment;
    catch err
        info.common_refs=[];
        info.common_refs_comment='';
    end 
    % 
    try
        info.common_gnd=parameters.CommonGround.NumericValue;
    info.common_gnd_comment=parameters.CommonGround.Comment;
    
    catch err
        info.common_gnd=[];
    info.common_gnd_comment='';
    
    end
    % 
    info.user_comment=parameters.UserComment.Value;
    % 
    info.audio_presentation=parameters.AudioSwitch.NumericValue;
    info.audio_presentation_comment=parameters.AudioSwitch.Comment;
    %
    info.filter_type=ecog_parameters.filter_type;
    info.bandpass_freq_low_stop=ecog_parameters.f_low_stop;
    info.bandpass_freq_low_pass=ecog_parameters.f_low_pass;
    info.bandpass_freq_high_stop=ecog_parameters.f_high_stop;
    info.bandpass_freq_high_pass=ecog_parameters.f_high_pass;
    info.bandpass_freq_range=[ecog_parameters.f_low_pass;ecog_parameters.f_high_pass]';
    % 
    info.noisy_channels=ecog_parameters.channels_noise;
    info.unselected_channels=ecog_parameters.channels_deselect;
    info.selected_channels=ecog_parameters.channels_selected;
    % 
    info.downsample_sampling_rate=ecog_parameters.samplingrate;
    % 
    if ~isempty(subject_op_info)
        info.session_info=session_info;
        info.operation_info=subject_op_info.op_info;
    else 
        info.session_info=[];
        info.operation_info=[];
    end 
    data=dat;
    if ~exist (strcat(save_path,subject_name))
        mkdir(strcat(save_path,subject_name));
    end 
    eval(strcat(subject_name,'_',session_name,'.data=dat')) ;
    eval(strcat(subject_name,'_',session_name,'.info=info'));
    %save(strcat(d(i).folder,'/',d(i).name),'data','info','-v7.3');
    save(strcat(save_path,subject_name,'/',subject_name,'_',session_name,'_crunched.mat'),strcat(subject_name,'_',session_name),'-v7.3');
    %clear states parameters 
    else 
       fprintf('corrupt experiment file \n'); 
    end
    %save(strcat(save_path,subject_name,'_',session_name,'_crunched.mat'),strcat(subject_name,'_',session_name),'-v7.3');
    %clear states parameters 
end


