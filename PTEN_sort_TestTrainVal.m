%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    seperate datasets for Train/Test/Val                   %
%                                                           %
%  this function balances into 70/20/10 by TMA              %
%     the patient level for each dataset seperately         %
%                                                           %
%       ..\PTEN                                             %
%             ..\5x                                         %  
%             ..\10x                                        %
%             ..\20x                                        %
%             ..\TMA                                        %
%                ..\train_allimgs                           %
%                ..\val_allimgs                             %
%                   ..\pos                                  % 
%                   ..\neg                                  %
%                ..\test_bypatient                          %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptDir = '/run/user/1000/gvfs/smb-share:server=nciis-p101.nci.nih.gov,share=group09/MIP/MIP/Stephanie Harmon/Queens_PTEN/classification/TMA';
saveDir = '/data/Stephanie_Harmon/PTEN';
%ptDir = [mainDir filesep 'raw' filesep 'decon'];
%saveDir = [mainDir filesep 'data' filesep 'decon'];

All_list = dir(ptDir);
All_list = All_list([All_list.isdir] & ~strcmpi({All_list.name},'.') & ~ strcmpi({All_list.name},'..'));
All_patients =  {All_list.name}';
exclude_list = {'TMA4_5_16_Complete';'TMA8_1_11_Partial';'TMA9_1_2_Complete'};
Others_patients = setdiff(All_patients,exclude_list);

% split to 10% test (+3 cases still needing processing)
% remaining 70% train and 20% val
rng(5); %set seed
Others_test  = randsample(Others_patients,round(numel(All_patients)*0.1)-3,false);
Others_remain = setdiff(Others_patients,Others_test);
Others_train = randsample(Others_remain,round(numel(All_patients)*0.7),false);
Others_val = setdiff(Others_remain, Others_train);
Others_test = cat(1,Others_test,exclude_list);
% parse over for each resolution
% create directories
for resi = 1:4
   if(resi == 1)
       res = '5x';
   elseif(resi == 2)
       res = '10x';
   elseif(resi == 3)
       res = '20x';
   elseif(resi==4)
       res = 'TMA';
   end
   disp(res)
       
        if(~exist([saveDir filesep res])) mkdir([saveDir filesep res]); end
        if(~exist([saveDir filesep res filesep 'test_bypatient']))  mkdir([saveDir filesep res filesep 'test_bypatient']); end
%         if(~exist([saveDir filesep res filesep 'train_allimgs']))   mkdir([saveDir filesep res filesep 'train_allimgs']); end
%         if(~exist([saveDir filesep res filesep 'train_obsv']))   mkdir([saveDir filesep res filesep 'train_obsv']); end
%         %if(~exist([saveDir filesep res filesep 'train_bypatient'])) mkdir([saveDir filesep res filesep 'train_bypatient']); end
%         if(~exist([saveDir filesep res filesep 'val_allimgs']))     mkdir([saveDir filesep res filesep 'val_allimgs']); end
%         %if(~exist([saveDir filesep res filesep 'val_bypatient']))   mkdir([saveDir filesep res filesep 'val_bypatient']); end
%         if(~exist([saveDir filesep res filesep 'val_allimgs' filesep 'pos']))     mkdir([saveDir filesep res filesep 'val_allimgs' filesep 'pos']); end
%         if(~exist([saveDir filesep res filesep 'val_allimgs' filesep 'neg']))     mkdir([saveDir filesep res filesep 'val_allimgs' filesep 'neg']); end
%         if(~exist([saveDir filesep res filesep 'train_obsv' filesep 'pos']))     mkdir([saveDir filesep res filesep 'train_obsv' filesep 'pos']); end
%         if(~exist([saveDir filesep res filesep 'train_obsv' filesep 'neg']))     mkdir([saveDir filesep res filesep 'train_obsv' filesep 'neg']); end

        test_pts = Others_test;
        train_pts = Others_train;
        val_pts = Others_val;

        %testing by patient
        disp('         ... starting send to testing')
        for patienti = 1:size(test_pts,1)
            test_files = dir([ptDir filesep test_pts{patienti} filesep res filesep '*.png']);
            if(contains(test_files(1).name,'Intact.png'))
                copyfile([ptDir filesep test_pts{patienti} filesep res],[saveDir filesep res filesep 'test_bypatient' filesep test_pts{patienti} '_pos'])
            elseif(~contains(test_files(1).name,'Intact.png'))
                copyfile([ptDir filesep test_pts{patienti} filesep res],[saveDir filesep res filesep 'test_bypatient' filesep test_pts{patienti} '_neg'])                
            end
        end

%         %training by patient
%         disp('         ... starting send to training')
%         if(size(train_pts,1)>0)
%             for patienti = 1:size(train_pts,1)
%                 %training by file
%                 train_files = dir([ptDir filesep train_pts{patienti} filesep res filesep '*.png']);
%                 
%                 for filei = 1:size(train_files,1)
%                     copyfile([train_files(filei).folder filesep train_files(filei).name],[saveDir filesep res filesep 'train_allimgs' filesep train_files(filei).name])
%                     
%                     if(contains(train_files(filei).name,'Intact.png'))
%                         copyfile([train_files(filei).folder filesep train_files(filei).name],[saveDir filesep res filesep 'train_obsv' filesep 'pos' filesep train_files(filei).name]);
%                     elseif(~contains(train_files(filei).name,'Intact.png'))
%                         copyfile([train_files(filei).folder filesep train_files(filei).name],[saveDir filesep res filesep 'train_obsv' filesep 'neg' filesep train_files(filei).name]);
%                     end
%                 end
%                 
%             end
%         end
% 
%         %validation by patient
%         disp('         ... starting send to validation')
%         if(size(val_pts,1)>0)
%             for patienti = 1:size(val_pts,1)
%                 %validation by file
%                 val_files = dir([ptDir filesep val_pts{patienti} filesep res filesep '*.png']);
%                 for filei = 1:size(val_files,1)
%                     if(contains(val_files(filei).name,'Intact.png'))
%                         copyfile([val_files(filei).folder filesep val_files(filei).name],[saveDir filesep res filesep 'val_allimgs' filesep 'pos' filesep val_files(filei).name]);
%                     elseif(~contains(val_files(filei).name,'Intact.png'))
%                         copyfile([val_files(filei).folder filesep val_files(filei).name],[saveDir filesep res filesep 'val_allimgs' filesep 'neg' filesep val_files(filei).name]);
%                     end
%                 end
%             end
%         end
% %    end
end

