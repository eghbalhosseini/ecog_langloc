% create subject info for langloc experiment 
save_path='/Users/eghbalhosseini/MyData/ecog-langloc/sub_operation_info/';
%% AMC 82
op_info=struct;
op_info.subject_name='AMC082';

op_info.handedness='R';
op_info.Language='N/A';
op_info.Age=51;
op_info.Sex='M';
op_info.GND=[127];
op_info.Ref=[128];
op_info.num_channels=128;
op_info.ecog_channels=1:126;
op_info.skull_eeg_channels=127:128;
op_info.microphone_channels=[192,129];
op_info.EMG_channels=[129:134];
op_info.visual_trigger='DigIO_1';
op_info.button_trigger='DigIO_2';
op_info.buzzer_trigger='DigIO_3';
op_info.audio_trigger='DigIO_4';
op_info.bad_channels = [25,26,13,101];
op_info.seizure_channels=[121,122,105,106,87,95,96,25,80];
labels= [arrayfun(@(x) sprintf('LG_%d',x),[1:64]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LFM_%d',x),[1:16]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LFS_%d',x),[1:16]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LFA_%d',x),[1:4]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LIP_%d',x),[1:12]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LIA_%d',x),[1:10]','uniformoutput',false);...
        arrayfun(@(x) sprintf('RIH_%d',x),[1:4]','uniformoutput',false);...
        'GND';'GND'];
op_info.channel_labels=labels;
% exp specific information
op.info.exp_name='MITLangloc';
op_info.task_file_location='/DAY3/MITLangloc/ECOG001/';
% 
op_info.time_of_runs={'18:56';'19:01';'19:10';'19:19';'19:29'};
op_info.run_descriptions={'MIT_Langloc_speedtest_1 - Aborted'; ... 
                      'MIT_Langloc_speedtest_1';...
                      'MIT_Langloc_slow_c1_1';...
                      'MIT_Langloc_slow_c1_2';...
                      'MIT_Langloc_slow_c1_3'};
                  
op_info.good_run=[0;1;1;1;1];

op_info.experimenter_notes={'Aborted due to issue with keyboard (patient inadvertently pressing other keys)'; '';'';'';''};

op_info.channel_labels=labels;
op_info.seizure_channels=[121,122,105,106,87,95,96];
op_info.exp_name='';
eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');


