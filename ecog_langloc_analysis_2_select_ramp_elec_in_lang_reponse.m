% selecting langauge responsive electrodes and add them to the subject
% info. 
% tested for subject 1 and validated with Terry's data . 
%% step 0: prepare the data 
clear all 
close all 
home 
%% 
data_path='/Users/eghbalhosseini/MyData/ecog-langloc/crunched/';
save_path='/Users/eghbalhosseini/MyData/ecog-langloc/crunched/';
analysis_path='/Users/eghbalhosseini/MyData/ecog-langloc/analysis/langloc_analysis_2_plot_chan_response';

d= dir([data_path,'/**/AMC097*_crunched.mat']);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
p_threshold=0.05;
do_print=1
%% 
sentences_electrode_with_langauge_accross_sessions=[];
jabberwocky_electrode_with_langauge_accross_sessions=[];
for i=2:length(d)
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
    sentence_gamma_band_ave_envelope=cellfun(@(x) x(1:12,gamma_band_index),{sentences.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    %sentence_gamma_band_ave_envelope=cellfun(@(x) x(1:12),{sentences.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*trial(column) structure
    sentence_gamma_band_ave_envelope=[sentence_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_ave_envelope,[3,2,1]));
    %append to langauge_channel*trial*words positions
    electrodes_with_language_response=sentence_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    sentences_electrode_with_langauge_accross_sessions=cat(2,sentences_electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    % 
    nonword_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'Jabberwocky'),info.word_type,'UniformOutput',false));
    nonwords=[data{nonword_trial_index}];
    nonword_gamma_band_ave_envelope=cellfun(@(x) x(1:12,gamma_band_index),{nonwords.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    %nonword_gamma_band_ave_envelope=cellfun(@(x) x(1:12),{nonwords.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    nonword_gamma_band_ave_envelope=[nonword_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    nonwords_gamma_band_ave_envelope_tensor=cell2mat(permute(nonword_gamma_band_ave_envelope,[3,2,1]));
    electrodes_response=nonwords_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    jabberwocky_electrode_with_langauge_accross_sessions=cat(2,jabberwocky_electrode_with_langauge_accross_sessions,electrodes_response);
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
end   

%% 
%% find ramping electrodes 
electrode_num=find(language_electrode);
sentence_all=[];
for i=1:length(electrode_num)
    channel_response=double(squeeze(sentences_electrode_with_langauge_accross_sessions(i,:,:)));
    word_position=repmat(1:size(channel_response,2),[size(channel_response,1),1]);
    [r_sentence,p_sent]=corr(mean(channel_response,1)',mean(word_position,1)','type','Spearman');
    sentence_all=[sentence_all;[r_sentence,p_sent]];
end 
ramp_electrodes=electrode_num(sentence_all(:,1)>0 & sentence_all(:,2)<0.01);
ramp_electrodes_location=(sentence_all(:,1)>0 & sentence_all(:,2)<0.01);
% 
electrode_with_ramp_across_sessions=sentences_electrode_with_langauge_accross_sessions(ramp_electrodes_location,:,:);
channel_ramp=zeros(size(language_electrode,1),1);
channel_ramp(ramp_electrodes)=1;
%% add ramping neurons 
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
    save(strcat(d(k).folder,'/',d(k).name),strcat(subject_name,'_',session_name),'-v7.3');
end
%% plot electrodes
close all;
electrode_num=find(language_electrode);
colors=cbrewer('div','RdYlBu',10);
num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;
for i=1:length(electrode_num)
    channel_response=double(squeeze(sentences_electrode_with_langauge_accross_sessions(i,:,:)));
    word_position=repmat(1:size(channel_response,2),[size(channel_response,1),1]);
    perturbed_word_position=word_position+.2*rand(size(word_position,1),size(word_position,2))-.1;
  
    fig1=figure(fix((i-1)/total_plots)+1);
    set(gcf,'Position',[-882 449 877 1325]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,electrode_num(i));
    a=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    % ah=violinPlot(channel_response, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0, ...
    %'color',  [1,.5,.5],'histOpt',1);
    hold on 
     h=scatter(perturbed_word_position(:)+.12,channel_response(:),7,'MarkerEdgeColor','none',...
              'MarkerFaceColor',colors(2,:))
    set(h,'MarkerFaceAlpha',.5)
    hAnnotation = get(h,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    
    jabberwocky_channel_response=double(squeeze(jabberwocky_electrode_with_langauge_accross_sessions(i,:,:)));
   %  ah1=violinPlot(jabberwocky_channel_response, 'histOri', 'left', 'widthDiv', [2 .9], 'showMM', 0, ...
    %'color',  [.7,.7,1],'histOpt',1);
 
    h=scatter(perturbed_word_position(:)-.12,jabberwocky_channel_response(:),7,'MarkerEdgeColor','none',...
              'MarkerFaceColor',colors(end-2,:))
    
    set(h,'MarkerFaceAlpha',.5)
    hAnnotation = get(h,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    e=plot(mean(word_position,1),mean(channel_response,1),...
        '-s','MarkerSize',1,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k','displayname','Sentences');
    e.Color=colors(2,:);
        

    f=plot(mean(word_position,1),mean(jabberwocky_channel_response,1),...
        '-s','MarkerSize',1,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k','displayname','Jabberwocky');
    
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:12])
    
    f.Color=colors(end-1,:);
    title(sub_title);
    set(a,'color','none','layer','top')
    new_bounds=quantile([channel_response(:);jabberwocky_channel_response(:)],99);
    
    a.YLim=[new_bounds(1),new_bounds(end)];
    if ~mod(i,total_plots) | i==length(electrode_num)
        legend('show','position',[.9,.1,.07,.05])
    if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path,'/',info.subject))
    end
    if do_print==1
        set(fig1,'PaperPosition',[.25 .25 8 6])
        print(fig1,'-bestfit', '-painters','-dpdf', strcat(analysis_path,'/',info.subject,'/',...
            info.subject,'_channel_activity','_fig_',num2str(fix((i-1)/num_rows)+2),'.pdf'));
    end
    end 
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
    sentence_gamma_band_envelope=cellfun(@(x) x(1:12,gamma_band_index),{sentences.signal_envelope_dowsample_parsed},'UniformOutput',false);
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
% sometimes, axes are plotted in a dark grey thats not exactly black (which
% % I find annoying). Make sure this doesnt happen.
% figHandles = get(groot, 'Children');
% for b= 1:length(figHandles)
%     axes = findobj(figHandles(b), 'type', 'axes');
%     for a = 1:length(axes),
%
%         if axes(a).YColor &amp;amp;amp;amp;lt; [1 1 1],
%             axes(a).YColor = [0 0 0];
%         end
%         if axes(a).XColor &amp;amp;amp;amp;lt; [1 1 1],
%             axes(a).XColor = [0 0 0];
%         end
%     end
%
%     print(figHandles(b), '-dpdf', strcat(subj_id,'_ramp_effect','.pdf'));
% end

% when you plotted several subplots but want them to have shared axes, use
% suplabel

% save to pdf
% see also export_fig from the file exchange

