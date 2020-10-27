% selecting langauge responsive electrodes and add them to the subject
% info. 
% tested for subject 1 and validated with Terry's data . 
%% step 0: prepare the data 
clear all 
%close all 
home 
%% 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-langloc/crunched/AMC092';
save_path='/Users/eghbalhosseiniasl1/MyData/ecog-langloc/crunched/AMC092';
d= dir([data_path,'/**/AMC092*_crunched.mat']);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
p_threshold=0.01;
%% 
electrode_with_langauge_accross_sessions=[];
electrode_nonwords_with_langauge_accross_sessions=[]
for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    try
    language_electrode=info.language_responsive_electrodes;
    end 
% step 1: extract electrodes with siginificant language response
    
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence 
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_ave_envelope=cellfun(@(x) x(1:12),{sentences.signal_ave_hilbert_zs_downsample_parsed},'UniformOutput',false);
    %sentence_gamma_band_ave_envelope=cellfun(@(x) x(1:12,gamma_band_index),{sentences.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*trial(column) structure
    sentence_gamma_band_ave_envelope=[sentence_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_ave_envelope,[3,2,1]));
    %append to langauge_channel*trial*words positions
    electrodes_with_language_response=sentence_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    electrode_with_langauge_accross_sessions=cat(2,electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    % 
    %  nonword trials 
    nonword_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'Jabberwocky'),info.word_type,'UniformOutput',false));
    nonwords=[data{nonword_trial_index}];
    nonword_gamma_band_ave_envelope=cellfun(@(x) x(1:12),{nonwords.signal_ave_hilbert_zs_downsample_parsed},'UniformOutput',false);
    %nonword_gamma_band_ave_envelope=cellfun(@(x) x(1:12,gamma_band_index),{nonwords.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    nonword_gamma_band_ave_envelope=[nonword_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    nonwords_gamma_band_ave_envelope_tensor=cell2mat(permute(nonword_gamma_band_ave_envelope,[3,2,1]));
    %words=~cell2mat(cellfun(@isempty,cellfun(@(x) strfind(x,'word'),[sentences.stimuli_type],'UniformOutput',false),'UniformOutput',false));
    nonwords_gamma_ave=nanmean(nonwords_gamma_band_ave_envelope_tensor,3);
    % add the new ones to the list 
    electrodes_non_nonwords_with_language_response=nonwords_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    electrode_nonwords_with_langauge_accross_sessions=cat(2,electrode_nonwords_with_langauge_accross_sessions,electrodes_with_language_response);
    
    
end   

%% find ramping electrodes 
electrode_num=find(language_electrode);
sentence_all=[];
for i=1:length(electrode_num)
    channel_response=double(squeeze(electrode_with_langauge_accross_sessions(i,:,:)));
    word_position=repmat(1:size(channel_response,2),[size(channel_response,1),1]);
    [r_sentence,p_sent]=corr(mean(channel_response,1)',mean(word_position,1)','type','Spearman');
    sentence_all=[sentence_all;[r_sentence,p_sent]];
end 
ramp_electrodes=electrode_num(sentence_all(:,1)>0 & sentence_all(:,2)<0.01);
ramp_electrodes_location=(sentence_all(:,1)>0 & sentence_all(:,2)<0.01);
% 
electrode_with_ramp_across_sessions=electrode_with_langauge_accross_sessions(ramp_electrodes_location,:,:);
channel_ramp=zeros(size(language_electrode,1),1);
channel_ramp(ramp_electrodes)=1;
%%  add the ramp electrode to back to the data 
for k=1:length(d)
    fprintf('adding ramp electrodes to %s \n', strcat(d(k).folder,'/',d(k).name));
    subj=load(strcat(d(k).folder,'/',d(k).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    info.ramp_electrodes=channel_ramp;
    subject_name=info.subject;
    session_name=info.session_name;
    
    eval(strcat(subject_name,'_',session_name,'.data=data;')) ;
    eval(strcat(subject_name,'_',session_name,'.info=info;'));
    %save(strcat(d(k).folder,'/',d(k).name),'data','info','-v7.3');
    save(strcat(save_path,subject_name,'_',session_name,'_crunched.mat'),strcat(subject_name,'_',session_name),'-v7.3');
end

%% plot electrodes
close all;

colors=cbrewer('qual','Set1',10);
num_rows=2;
num_columns=1;
total_plots=num_rows*num_columns;

for i=1:length(ramp_electrodes)
    electrode_response=squeeze(electrode_with_ramp_across_sessions(i,:,:));
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,ramp_electrodes(i));
    word_position=repmat(1:size(electrode_response,2),[size(electrode_response,1),1]);
    perturbed_word_position=word_position+.2*rand(size(electrode_response,1),size(electrode_response,2))-.1;
   
  
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[-1685 519 492 1131]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,ramp_electrodes(i));
    a=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
     %ah=violinPlot(channel_response, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0, ...
    %'color',  mat2cell(colors(1, : ), 1),'histOpt',0);
    hold on 
    h=scatter(perturbed_word_position(:),channel_response(:),4);
    set(h,'CData',colors(6,:),'MarkerFacecolor',colors(1,:),'MarkerFaceAlpha',.3)
   
    e=plot(mean(word_position,1),mean(channel_response,1),...
        '-s','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k');
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:size(word_position,2)])
    xlabel('word position')
    ylabel('high \gamma z-score')
    e.Color=[0,0,0];
    title(sub_title);
    set(a,'color','none','layer','top')
    
end  



%% plot ramp electrodes 
figure;
set(gcf,'OuterPosition',[-1615 441 694 1337]);
plot_width=.9/length(d);
plot_length=.87;
for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
    
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:12),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_envelope_tensor=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    if i==1
        a=mean(sentence_gamma_band_envelope_tensor,3);
        b=mean(a,2);
        [~,idx]=sort(b);
        axes('position',[.05,0.05+plot_width*(i-1),plot_length,plot_width]);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        a=mean(sentence_gamma_band_envelope_tensor,3);
        imagesc(a(idx,:));
        set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
        handles = colorbar;
        handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'z-score';
        drawnow;
        xlabel('word position'); ylabel('Electrode');
    else 
    axes('position',[.05,0.05+plot_width*(i-1)+0.005*(i-1),plot_length,plot_width])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    a=mean(sentence_gamma_band_envelope_tensor,3);
    imagesc(a(idx,:));
    set(gca, 'ydir', 'normal','box','off','xtick',[],'xticklabel','');    
    handles = colorbar;
    handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
    end 
    
end
print(gcf, '-depsc', strcat(info.subject,'_sentence_ramp_elec','.eps'));
%% 
figure;
set(gcf,'OuterPosition',[-1615 441 694 1337]);
plot_width=.9/length(d);
plot_length=.87;
for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
    
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_envelope_tensor=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    if i==1
        a=mean(sentence_gamma_band_envelope_tensor,3);
        b=mean(a,2);
        [~,idx]=sort(b);
        axes('position',[.05,0.05+plot_width*(i-1),plot_length,plot_width]);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        a=mean(sentence_gamma_band_envelope_tensor,3);
        imagesc(a(idx,:));
        set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
        handles = colorbar;
        handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'z-score';
        drawnow;
        xlabel('word position'); ylabel('Electrode');
    else 
    axes('position',[.05,0.05+plot_width*(i-1)+0.005*(i-1),plot_length,plot_width])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    a=mean(sentence_gamma_band_envelope_tensor,3);
    imagesc(a(idx,:));
    set(gca, 'ydir', 'normal','box','off','xtick',[],'xticklabel','');    
    handles = colorbar;
    handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
    end 
    
end
print(gcf, '-depsc', strcat(info.subject,'_Jaberwocky_res_ramp','.eps'));
%% 
%% 
figure;
set(gcf,'OuterPosition',[-1615 441 694 1337]);
plot_width=.9/length(d);
plot_length=.87;
for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
    
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_envelope_tensor=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    if i==1
        a=mean(sentence_gamma_band_envelope_tensor,3);
        b=mean(a,2);
        [~,idx]=sort(b);
        axes('position',[.05,0.05+plot_width*(i-1),plot_length,plot_width]);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        a=mean(sentence_gamma_band_envelope_tensor,3);
        imagesc(a(idx,:));
        set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
        handles = colorbar;
        handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'z-score';
        drawnow;
        xlabel('word position'); ylabel('Electrode');
    else 
    axes('position',[.05,0.05+plot_width*(i-1)+0.005*(i-1),plot_length,plot_width])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    a=mean(sentence_gamma_band_envelope_tensor,3);
    imagesc(a(idx,:));
    set(gca, 'ydir', 'normal','box','off','xtick',[],'xticklabel','');    
    handles = colorbar;
    handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
    end 
    
end
print(gcf, '-depsc', strcat(info.subject,'_words_res_ramp','.eps'));
%% 
figure;
set(gcf,'OuterPosition',[-1615 441 694 1337]);
plot_width=.9/length(d);
plot_length=.87;
for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
    
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_envelope_tensor=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    if i==1
        a=mean(sentence_gamma_band_envelope_tensor,3);
        b=mean(a,2);
        [~,idx]=sort(b);
        axes('position',[.05,0.05+plot_width*(i-1),plot_length,plot_width]);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        a=mean(sentence_gamma_band_envelope_tensor,3);
        imagesc(a(idx,:));
        set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
        handles = colorbar;
        handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'z-score';
        drawnow;
        xlabel('word position'); ylabel('Electrode');
    else 
    axes('position',[.05,0.05+plot_width*(i-1)+0.005*(i-1),plot_length,plot_width])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    a=mean(sentence_gamma_band_envelope_tensor,3);
    imagesc(a(idx,:));
    set(gca, 'ydir', 'normal','box','off','xtick',[],'xticklabel','');    
    handles = colorbar;
    handles.Position=[0.93 0.05+plot_width*(i-1)+0.005*(i-1) .02 0.09];
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
    end 
    
end
print(gcf, '-depsc', strcat(info.subject,'_nonwords_res_ramp','.eps'));


%% 
%% plot mean language electrodes
session_sentence_gamma_band_envelope_tensor=[];
session_words_gamma_band_envelope_tensor=[];
session_nonwords_gamma_band_envelope_tensor=[];
session_jabberwocky_gamma_band_envelope_tensor=[];
for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_envelope_tensor=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    session_sentence_gamma_band_envelope_tensor=cat(3,session_sentence_gamma_band_envelope_tensor,sentence_gamma_band_envelope_tensor);
    %
        % step 1: extract electrodes with siginificant language response
    nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    % nonwords
    nonwordss=[data{nonwords_trial_index}];
    nonwords_gamma_band_envelope=cellfun(@(x) x(1:8),{nonwordss.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    nonwords_gamma_band_envelope=[nonwords_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    nonwords_gamma_band_envelope_tensor=cell2mat(permute(nonwords_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    nonwords_gamma_band_envelope_tensor=nonwords_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    session_nonwords_gamma_band_envelope_tensor=cat(3,session_nonwords_gamma_band_envelope_tensor,nonwords_gamma_band_envelope_tensor);
    % 
    % step 1: extract electrodes with siginificant language response
    jabberwocky_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
    % jabberwocky
    jabberwockys=[data{jabberwocky_trial_index}];
    jabberwocky_gamma_band_envelope=cellfun(@(x) x(1:8),{jabberwockys.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    jabberwocky_gamma_band_envelope=[jabberwocky_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    jabberwocky_gamma_band_envelope_tensor=cell2mat(permute(jabberwocky_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    jabberwocky_gamma_band_envelope_tensor=jabberwocky_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    session_jabberwocky_gamma_band_envelope_tensor=cat(3,session_jabberwocky_gamma_band_envelope_tensor,jabberwocky_gamma_band_envelope_tensor);
    % 
    % step 1: extract electrodes with siginificant language response
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    % words
    wordss=[data{words_trial_index}];
    words_gamma_band_envelope=cellfun(@(x) x(1:8),{wordss.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    words_gamma_band_envelope=[words_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    words_gamma_band_envelope_tensor=cell2mat(permute(words_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    words_gamma_band_envelope_tensor=words_gamma_band_envelope_tensor(find(language_electrode),:,:);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    session_words_gamma_band_envelope_tensor=cat(3,session_words_gamma_band_envelope_tensor,words_gamma_band_envelope_tensor);
end
%% 
temp=[mean(session_sentence_gamma_band_envelope_tensor,3),...
    mean(session_jabberwocky_gamma_band_envelope_tensor,3),...
    mean(session_nonwords_gamma_band_envelope_tensor,3),...
    mean(session_words_gamma_band_envelope_tensor,3)];
min_ax=min(temp(:));
max_ax=max(temp(:));
figure 
plot_width=.18;
plot_length=.8;
set(gcf,'OuterPosition',[-1615 441 694 1337]);
a=mean(session_sentence_gamma_band_envelope_tensor,3);
b=mean(a,2);
[~,idx]=sort(b);
axes('position',[.08,.05,plot_length,plot_width]);
title('sentences')
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(a(idx,:),[min_ax,max_ax]);
set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
title('sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Position=[0.9 0.05 .02 plot_width];
handles.Label.String = 'z-score';
drawnow;
xlabel('word position'); ylabel('Electrode');

a=mean(session_nonwords_gamma_band_envelope_tensor,3);
axes('position',[.08,.28,plot_length,plot_width]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(a(idx,:),[min_ax,max_ax]);
set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
title('nonwords');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Position=[0.9 0.28 .02 plot_width];
handles.Label.String = 'z-score';
drawnow;
xlabel(''); ylabel('');

a=mean(session_jabberwocky_gamma_band_envelope_tensor,3);
axes('position',[.08,.5,plot_length,plot_width]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(a(idx,:),[min_ax,max_ax]);
set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
title('jabberwocky');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Position=[0.9 0.5 .02 plot_width];
handles.Label.String = 'z-score';
drawnow;
xlabel(''); ylabel('');


% 

a=mean(session_words_gamma_band_envelope_tensor,3);
axes('position',[.08,.75,plot_length,plot_width]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(a(idx,:),[min_ax,max_ax]);
set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2),'xticklabel',[0:8]);
title('words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Position=[0.9 0.75 .02 plot_width];
handles.Label.String = 'z-score';
drawnow;
xlabel(''); ylabel('');
print(gcf, '-depsc', strcat(info.subject,'_ramp_condition_comparison','.eps'));
%% 
S=double(mean(session_sentence_gamma_band_envelope_tensor,3));
W=double(mean(session_words_gamma_band_envelope_tensor,3));
N=double(mean(session_nonwords_gamma_band_envelope_tensor,3));
J=double(mean(session_jabberwocky_gamma_band_envelope_tensor,3));
cval={'S' 'W' 'J','N'};
cind=[ones(size(S,1),1);2*ones(size(W,1),1);3*ones(size(J,1),1);4*ones(size(N,1),1)];
c=cval(cind);
x=1:size(S,2);
y=[S;W;J;N];

clear g
g=gramm('x',x,'y',y,'color',c);
g.set_color_options('map','d3_20');
g.geom_hline();
g.set_order_options('color',[1,2,3,4]);
g.stat_summary();
g.set_title('mean activity across all channels ');
g.axe_property('Xtick',0.0:135:size(sentence_gamma_band_envelope_tensor,2));
g.axe_property('Xticklabel',[0:8]);
g.set_names('x','word position','y','z-score','color','stimuli');
figure('Position',[100 100 800 550]);
g.draw();
g.export('file_name',strcat(info.subject,'_ramp_condition_comparison_ave'),'file_type','pdf');
%print(gcf, '-depsc', strcat(info.subject,'_condition_comparison_ave','.eps'));





