% selecting langauge responsive electrodes and add them to the subject
% info. 
% tested for subject 1 and validated with Terry's data . 
%% step 0: prepare the data 
clear all 
close all 

data_path='/Users/eghbalhosseini/MyData/ecog-langloc/';
save_path='/Users/eghbalhosseini/MyData/ecog-langloc/crunched/';
d= dir([data_path,'/**/AMC097*_crunched.mat']);
subject_id='AMC082';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);

fprintf(' %d .mat files were found \n', length(d));
num_of_permutation=1000; 
p_threshold=0.05;
%% 
trial_condition=[];
trial_hilbert_ave=[];
trial_gamma_ave=[];
valid_chann=[];
%% 

for k=1:length(d)
    fprintf('adding %s from %s \n',d(k).name, strcat(d(k).folder,'/',d(k).name));
    subj=load(strcat(d(k).folder,'/',d(k).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    valid_chann=[valid_chann,info.valid_channels];
    
    % step 1: compute mean across the word positions in each trial 

    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'Sentences'),info.word_type,'UniformOutput',false));
    % sentence 
    sentences=[data{sentence_trial_index}];
    sentence_hilbert_band_ave_envelope=cellfun(@(x) x(1:12),{sentences.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    sentence_hilbert_band_ave_envelope=[sentence_hilbert_band_ave_envelope{:,:}];
    sentence_hilbert_band_ave_envelope_tensor=cell2mat(permute(sentence_hilbert_band_ave_envelope,[3,2,1]));
    %  nonword trials 
    nonword_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'Jabberwocky'),info.word_type,'UniformOutput',false));
    nonwords=[data{nonword_trial_index}];
    nonword_hilbert_band_ave_envelope=cellfun(@(x) x(1:8),{nonwords.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    nonword_hilbert_band_ave_envelope=[nonword_hilbert_band_ave_envelope{:,:}];
    nonwords_hilbet_band_ave_envelope_tensor=cell2mat(permute(nonword_hilbert_band_ave_envelope,[3,2,1]));
    % add the new ones to the list 
    trial_hilbert_ave=[trial_hilbert_ave,[nanmean(sentence_hilbert_band_ave_envelope_tensor,3),nanmean(nonwords_hilbet_band_ave_envelope_tensor,3)]];
    trial_condition=[trial_condition,[nanmean(sentence_hilbert_band_ave_envelope_tensor,3)*0+1,nanmean(nonwords_hilbet_band_ave_envelope_tensor,3)*0-1]];
    fprintf('added %s from %s \n',d(k).name, strcat(d(k).folder,'/',d(k).name));
end
    
%% step 2: compute a correlation between trial means and vector of condition labels )
% sentences = 1, nonword-lists =-1 
num_repeats=1000;
[RHO_hilbert,~] = corr(double(transpose(trial_hilbert_ave)),double(transpose(trial_condition)),'Type','Spearman','rows','complete');

rho_hilbert_original=diag(RHO_hilbert);
rho_hilbert_positive=rho_hilbert_original>0;
% step 3: random permutation of conditions labels, repeat 1000 times 
rho_hilbert_permuted=[];
fprintf('permutation :');
for k=1:num_repeats
    if ~mod(k,100), fprintf(' %d ',k ); end
    random_index=randperm(size(trial_condition,2));
    [RHO_hilbert_rand,~] = corr(double(transpose(trial_hilbert_ave)),double(transpose(trial_condition(:,random_index))),'Type','Spearman','rows','complete'); 
    rho_hilbert_permuted=[rho_hilbert_permuted,diag(RHO_hilbert_rand)];
end

%% step 4 : compute fraction of correlations in step 3 that produce higher correlation that step 2. 
% step 4 : compute fraction of correlations in step 3 that produce higher correlation that step 2. 
p_fraction=sum(rho_hilbert_permuted>repmat(rho_hilbert_original,[1,num_repeats]),2)./size(rho_hilbert_permuted,2);
p_significant=p_fraction<p_threshold;
channel_hilbert_significant=p_significant.*rho_hilbert_positive;
if size(info.valid_channels,1)<size(info.valid_channels,2)
    info.valid_channels=transpose(info.valid_channels);
end 
channel_hilbert_significant=channel_hilbert_significant.*info.valid_channels;
%channel_hilbert_significant(info.unselected_channels)=0;
fprintf('\n hilbert language electrodes: ');
fprintf(' %d  ',find(channel_hilbert_significant)' );
fprintf('\n');
%% step 5: add the language electrode to back to the data 
% step 5: add the language electrode to back to the data 
for pp=1:length(d)
    fprintf('adding language electrodes to %s \n', strcat(d(pp).folder,'/',d(pp).name));
    subj=load(strcat(d(pp).folder,'/',d(pp).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    info.language_responsive_electrodes_hilbert_odd=channel_hilbert_significant;
    if size(info.valid_channels,1)<size(info.valid_channels,2)
    info.valid_channels=transpose(info.valid_channels);
    end 
    subject_name=info.subject;
    session_name=info.session_name;
    eval(strcat(subject_name,'_',session_name,'.data=data;'));
    eval(strcat(subject_name,'_',session_name,'.info=info;'));
    save(strcat(save_path,subject_name,'_',session_name,'_crunched.mat'),strcat(subject_name,'_',session_name),'-v7.3');
end 
%% 
clear all 
close all 
data_path='/Users/eghbalhosseini/MyData/ecog-langloc';
save_path='/Users/eghbalhosseini/MyData/ecog-langloc/crunched/';
analysis_path='/Users/eghbalhosseini/MyData/ecog-langloc/analysis/analysis_0_plot_language_responsive_elec/';
subject_id='AMC082';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
% 
session_sentence_hilbert_band_envelope_pre_trial_tensor=[];
session_sentence_hilbert_band_envelope_tensor=[];
% 

% 
session_nonwords_hilbert_band_envelope_pre_trial_tensor=[];
session_nonwords_hilbert_band_envelope_tensor=[];
% 
for i=1:length(d)
    
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.language_responsive_electrodes_hilbert_odd;
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    %%%%%%%%%%%%%%%%% sentence
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    sentences=[data{sentence_trial_index}];
    %word_length=sentences(1).signal_range_downsample(1,2)-sentences(1).signal_range_downsample(1,1)+1;
    %
    hilbert_band_envelope=cellfun(@(x) x(1:16),{sentences.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];

    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2])); 
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]); 
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{sentences.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial_tensor=[hilbert_band_envelope_pre_trial{:,:}];
    session_sentence_hilbert_band_envelope_pre_trial_tensor=cat(2,session_sentence_hilbert_band_envelope_pre_trial_tensor,hilbert_band_envelope_pre_trial_tensor(find(language_electrode),:,:));
    session_sentence_hilbert_band_envelope_tensor=cat(2,session_sentence_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor(find(language_electrode),:,:));
    %%%%%%%%%%%%%%%%% words
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    words=[data{words_trial_index}];
    %
    hilbert_band_envelope=cellfun(@(x) x(1:10),{words.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2])); 
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]); 
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{words.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial_tensor=[hilbert_band_envelope_pre_trial{:,:}];
    session_words_hilbert_band_envelope_pre_trial_tensor=cat(2,session_words_hilbert_band_envelope_pre_trial_tensor,hilbert_band_envelope_pre_trial_tensor(find(language_electrode),:,:));
    session_words_hilbert_band_envelope_tensor=cat(2,session_words_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor(find(language_electrode),:,:));
    %%%%%%%%%%%%%%%% nonwords
    nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    nonwords=[data{nonwords_trial_index}];
    %
    hilbert_band_envelope=cellfun(@(x) x(1:10),{nonwords.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2])); 
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]); 
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial_tensor=[hilbert_band_envelope_pre_trial{:,:}];
    session_nonwords_hilbert_band_envelope_pre_trial_tensor=cat(2,session_nonwords_hilbert_band_envelope_pre_trial_tensor,hilbert_band_envelope_pre_trial_tensor(find(language_electrode),:,:));
    session_nonwords_hilbert_band_envelope_tensor=cat(2,session_nonwords_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor(find(language_electrode),:,:));
    %%%%%%%%%%%%%%% jabberwocky
    jabberwocky_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
    jabberwocky=[data{jabberwocky_trial_index}];
    %
    hilbert_band_envelope=cellfun(@(x) x(1:10),{jabberwocky.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2])); 
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]); 
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial_tensor=[hilbert_band_envelope_pre_trial{:,:}];
    session_jabberwocky_hilbert_band_envelope_pre_trial_tensor=cat(2,session_jabberwocky_hilbert_band_envelope_pre_trial_tensor,hilbert_band_envelope_pre_trial_tensor(find(language_electrode),:,:));
    session_jabberwocky_hilbert_band_envelope_tensor=cat(2,session_jabberwocky_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor(find(language_electrode),:,:));
    clear subj data
    
end 
% 
close all 
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
colors=flipud(cbrewer('qual','Set1',5));
for i=1:length(language_electrode_num)
    % extract sentences 
    electrode_sentence_pre_trial_response=double(squeeze(session_sentence_hilbert_band_envelope_pre_trial_tensor(i,:,:)));
    electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_tensor(i,:,:)));
    electrode_nonwords_pre_trial_response=double(squeeze(session_nonwords_hilbert_band_envelope_pre_trial_tensor(i,:,:)));
    electrode_nonwords_response=double(squeeze(session_nonwords_hilbert_band_envelope_tensor(i,:,:)));
    electrode_words_pre_trial_response=double(squeeze(session_words_hilbert_band_envelope_pre_trial_tensor(i,:,:)));
    electrode_words_response=double(squeeze(session_words_hilbert_band_envelope_tensor(i,:,:)));
    electrode_jabber_pre_trial_response=double(squeeze(session_jabberwocky_hilbert_band_envelope_pre_trial_tensor(i,:,:)));
    electrode_jabber_response=double(squeeze(session_jabberwocky_hilbert_band_envelope_tensor(i,:,:)));
    % sentence
    % 
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[ 1442,-129,1068,1269])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    % 
    x=[-size(electrode_sentence_pre_trial_response,2)+1:0,1:size(electrode_sentence_response,2)];
    y=[nanmean(electrode_sentence_pre_trial_response,1),nanmean(electrode_sentence_response,1)];
    z=[nanstd(electrode_sentence_pre_trial_response,[],1),nanstd(electrode_sentence_response,[],1)]./sqrt(size(electrode_sentence_response,1));
    bl = boundedline(x, y, z, 'cmap', colors(1,:),'alpha');
    bl.LineWidth=2;
    hold on 
    bl.DisplayName='sentences';
    % 
    x=[-size(electrode_nonwords_pre_trial_response,2)+1:0,1:size(electrode_nonwords_response,2)];
    y=[nanmean(electrode_nonwords_pre_trial_response,1),nanmean(electrode_nonwords_response,1)];
    z=[nanstd(electrode_nonwords_pre_trial_response,[],1),nanstd(electrode_nonwords_response,[],1)]./sqrt(size(electrode_nonwords_response,1));
    bl = boundedline(x, y, z, 'cmap', colors(2,:),'alpha');
    bl.LineWidth=2;
    bl.DisplayName='Nonwords';
    % 
    x=[-size(electrode_words_pre_trial_response,2)+1:0,1:size(electrode_words_response,2)];
    y=[nanmean(electrode_words_pre_trial_response,1),nanmean(electrode_words_response,1)];
    z=[nanstd(electrode_words_pre_trial_response,[],1),nanstd(electrode_words_response,[],1)]./sqrt(size(electrode_words_response,1));
    bl = boundedline(x, y, z, 'cmap', colors(3,:),'alpha');
    bl.LineWidth=2;
    bl.DisplayName='words';
    % 
    x=[-size(electrode_jabber_pre_trial_response,2)+1:0,1:size(electrode_jabber_response,2)];
    y=[nanmean(electrode_jabber_pre_trial_response,1),nanmean(electrode_jabber_response,1)];
    z=[nanstd(electrode_jabber_pre_trial_response,[],1),nanstd(electrode_jabber_response,[],1)]./sqrt(size(electrode_jabber_response,1));
    bl = boundedline(x, y, z, 'cmap', colors(4,:),'alpha');
    bl.LineWidth=2;
    bl.DisplayName='jabberwocky';
    axis tight
    e1=arrayfun(@(x,y) plot([x,x],get(gca,'ylim'),'k--'),[0:word_length:8*word_length]);
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e1);hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
   
    % nonwords 
   
    % 
    ax.XAxis.Visible = 'on'; 
    set(ax,'ydir', 'normal','box','off');
    title(sub_title);
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        legend('show','position',[.92,.1,.07,.05])
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        
        print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_SWJN_response_pre_post',num2str(p),'_v2.eps')); 
    end 
    
end