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
    addpath(genpath('~/MyCodes/ecog-langloc'));
end 
    
%% 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-langloc';
save_path='/Users/eghbalhosseiniasl1/MyData/ecog-langloc/crunched/';
if ~exist(save_path)
    mkdir(save_path);
end
d= dir([data_path,'/**/ECOG*.dat']);
fprintf(' %d .dat files were found \n', length(d))
%% 
d_op_info= dir([data_path,'/**/*operation_info.mat']);
fprintf(' %d op_info files were found \n', length(d_op_info));

%% 
for i=2%:5%length(d)
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

    % 
    if ~isempty(subject_session_info)
        ecog_ch=subject_session_info.op_info.ecog_ch;
        ref_ch=subject_session_info.op_info.Ref;
        gnd_ch=subject_session_info.op_info.GND;
        bad_ch=subject_session_info.session_info.bad_channels{:};
        all_channels=ecog_ch;
        channels=-ecog_ch;
        %channels = setdiff(ecog_ch,unique([ref_ch,gnd_ch,bad_ch]));
        channels(setdiff(ecog_ch,unique([ref_ch,gnd_ch,bad_ch])))= -1*channels(setdiff(ecog_ch,unique([ref_ch,gnd_ch,bad_ch])));
    else
        ecog_ch=[]; gnd_ch=[]; bad_ch=[];
        ref_ch=parameters.RefChList.NumericValue;
        channels = setdiff(1:parameters.SourceCh.NumericValue',unique([ref_ch,gnd_ch,bad_ch]));
    end
    subject_data_dr=strcat(d(i).folder,'/',d(i).name);
    [~,~,parameters]=load_bcidat(subject_data_dr,[0,0]);
    cfg = [];
    cfg.dataset = subject_data_dr;
    cfg.trialfun = 'read_langloc_trials';
    cfg.trialdef.prestim = .5;             % pre-trial sample, in seconds
    cfg.trialdef.poststim = 0;            % post-trial sample, in seconds
    cfg.trialdef.eventvalue = 0;  % trial code (to identify the relevant time-points)
    cfg = ft_definetrial(cfg);
    cfg.channel=all_channels;
    cfg.padding='data';
    data=ft_preprocessing(cfg); 
    
    cfg=[];
    cfg.channel=all_channels;
      
    % high pass 
    disp('High-Pass filtering: ');
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;  % hpfilt_param.Wp;
    cfg.hpfilttype = 'but';      % default: butterworth
    cfg.hpfiltdir = 'twopass';  % default: forward and reverse filtering
    % cfg.detrend = 'yes';
    data = ft_preprocessing(cfg,data);
    cfg= [];
    cfg.channel=channels;
    cfg.viewmode = 'vertical';
    cfg = ft_databrowser(cfg, data);
    data = ft_rejectartifact(cfg, data);   
    % 

    disp('Common-average referencing:');
    channelDim=1;
    a=cellfun(@str2num,data.hdr.label,'UniformOutput',false);
    gAmp=1;
    if gAmp
    for trialInd = 1:length(data.trial)
          data2 = nan*zeros(size(data.trial{trialInd}));
            nAmps = ceil(max(data.hdr.nChans) / 16);  % number of separate amplifiers
            for ampInd = 1:nAmps
                firstChan  = (ampInd-1)*16+1;
                lastChan = ampInd*16;
                inds = find(and([a{:}]>=firstChan , [a{:}]<=lastChan));
                data2(inds,:) = ...
                    car(data.trial{trialInd}(inds,:),[],channelDim);              
            end
            data.trial{trialInd} = data2;
        end
    else
        for trialInd = 1:length(data.trial)
            data2 = car(data.trial{trialInd},[],channelDim);
            data.trial{trialInd} = data2;
        end
    end
    clear data2  
    disp(' ');
    
    disp('Notch filtering:' );
    cfg = [];
    cfg.channel=channels;
    cfg.dftfilter = 'yes';
    cfg.dftfreq = [60 120 180 240];
    data = ft_preprocessing(cfg,data);
    
    cfg= [];
    cfg.ylim=[-50,50];
    %cfg.colorgroups=1:(data.hdr.nChans);
    %colors=cbrewer('qual', 'Dark2', data.hdr.nChans);
    %cfg.channelcolormap=colors(randperm(size(colors,1)),:);
    cfg.linewidth=1.0;
    cfg.channel=channels;
    cfg.viewmode = 'vertical';
    cfg = ft_databrowser(cfg, data);
    data = ft_rejectartifact(cfg, data); 
    
    

    
    
    
   %% Bandpass filtering: Gerv's method  |  SANITY CHECKED %%
    bands = {'Alpha', [8 12];
    'Beta', [18 25];
    'LowGamma', [36 52];
    'HighGamma', [70 170]};
    if 1
        disp('Computing envelope: ');
        for bandInd = 1:size(bands,1)        
        %% Bandpass filter %%
        disp([bands{bandInd,1}, ' band:']);
        cfg = [];
        cfg.channel=channels;
        cfg.bpfilter = 'yes';
        cfg.bpfreq = bands{bandInd,2};      % high-pass and low-pass frequencies, from bfilt_param
        cfg.bpfilttype = 'but';             % default: IIR butterworth
        cfg.bpfiltdir = 'twopass';          % default: forward and reverse filtering
        data2 = ft_preprocessing(cfg,data);

        %% Hilbert transform: compute envelope (amplitude of analytic signal) %%
        cfg = [];
        cfg.hilbert = 'abs';
        env = ft_preprocessing(cfg,data2);       
        clear data2
        
        %% Smooth signal (low-pass filter) and downsample %%
        cfg = [];
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 100;           % lpfilt_param.Wp
        cfg.lpfilttype = 'but';     % default: butterworth
        cfg.lpfiltdir = 'twopass';  % default: forward and reverse filtering
        env2 = ft_preprocessing(cfg,env);
        clear env
        env = env2;
        clear env2
                
        cfg = [];
        cfg.resamplefs = 300;   % downsampling frequency
        cfg.detrend = 'no';             % this value has no default so it has to be set
        env2 = ft_resampledata(cfg,env);
        cfg= [];
        cfg.channel=channels;
        cfg.viewmode = 'vertical';
        cfg = ft_databrowser(cfg, env2);
        env2 = ft_rejectartifact(cfg, env2);
        end 
    end
    

    
    else 
        fprintf('bad dataset \n')
    end 

% 


end


