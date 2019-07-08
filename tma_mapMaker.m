srcl = 'M:/'; %

mainDir = [srcl 'Stephanie Harmon/Queens_PTEN/classification/TMA']; %'V:\TCGABLCA'
csvLoc = [srcl 'Stephanie Harmon/Queens_PTEN/predictions/TMA/csv'];

%val
mapLoc = [srcl 'Stephanie Harmon/Queens_PTEN/predictions/TMA/maps/images'];
mapLoc2 = [srcl 'Stephanie Harmon/Queens_PTEN/predictions/TMA/maps/matfiles'];
preds_5x  = readtable([csvLoc filesep 'val_results_by_img_5x_obs_5loss_07032019-1250.pkl_07052019-1216.csv'],'Delimiter',',');
preds_10x = readtable([csvLoc filesep 'val_results_by_img_10x_obs_5loss_07032019-1345.pkl_07052019-1111.csv'],'Delimiter',',');
preds_20x = readtable([csvLoc filesep 'val_results_by_img_20x_obs_one5loss_07052019-1019.pkl_07052019-1116.csv'],'Delimiter',',');


%find unique patients in 5x, 10x, 20x
list_5x = table2cell(preds_5x); list_5x = list_5x(:,2:5);
list_10x = table2cell(preds_10x); list_10x = list_10x(:,2:5);
list_20x = table2cell(preds_20x); list_20x = list_20x(:,2:5);


%5x
for i = 1:size(list_5x,1)
    filei = strsplit(list_5x{i,1},'_');
    list_5x{i,5} = strjoin(filei(1:3),'_');
    list_5x{i,6} = filei{4};
    list_5x{i,7} = strrep(filei{5},'.png','');
end
uniq_5x = unique(list_5x(:,5));
%10x
for i = 1:size(list_10x,1)
    filei = strsplit(list_10x{i,1},'_');
    list_10x{i,5} = strjoin(filei(1:3),'_');
    list_10x{i,6} = filei{4};
    list_10x{i,7} = filei{5};
end
uniq_10x = unique(list_10x(:,5));
%20x
for i = 1:size(list_20x,1)
    filei = strsplit(list_20x{i,1},'_');
    list_20x{i,5} = strjoin(filei(1:3),'_');
    list_20x{i,6} = filei{4};
    list_20x{i,7} = filei{5};
    list_20x{i,8} = filei{6};
end
uniq_20x = unique(list_20x(:,5));

uniq_all = cat(1,uniq_5x,uniq_10x,uniq_20x);
slides_unique = unique(uniq_all);

for casei = 1:size(slides_unique,1)
    disp(slides_unique{casei})
    slide_Data_5x  = list_5x(strcmpi(list_5x(:,5),slides_unique{casei}),:);
    slide_Data_10x = list_10x(strcmpi(list_10x(:,5),slides_unique{casei}),:);
    slide_Data_20x = list_20x(strcmpi(list_20x(:,5),slides_unique{casei}),:);
    tma_id = slides_unique{casei};
    tma_outcome = unique(slide_Data_5x(:,7)); tma_outcome = tma_outcome{1};
    
    %load entire image from TMA
    imgfind = dir([mainDir filesep tma_id '_' tma_outcome filesep 'TMA' filesep '*.png']);
    tma_img = imread([imgfind.folder filesep imgfind.name]);
    tma_img = imresize(tma_img,2);
    bw_img = double(rgb2gray(tma_img))./255;
    tissue = zeros(size(bw_img));
    tissue(find(bw_img<0.9)) = 1;
    tissue = imfill(tissue,'holes');
    se = strel('disk',25);
    tissue = imclose(tissue,se);
    tissue = bwareaopen(tissue, 10);
    tissue_sm = imresize(tissue,0.01,'method','nearest');
    
    imgsize = size(bw_img);
    %imgsize = ceil(imgsize_orig.*0.01);

    %roi_mask = zeros(imgsize_orig);
    count_map_5x = zeros(imgsize);
    count_map_10x = zeros(imgsize);
    count_map_20x = zeros(imgsize);
    count_map_40x = zeros(imgsize);

    out_map_5x = zeros(imgsize);
    out_map_10x = zeros(imgsize);
    out_map_20x = zeros(imgsize);
    out_map_40x = zeros(imgsize);

    all_mask = zeros(imgsize);
    


    %5x is 4 pixels
    %10x is 2 pixels
    %20x is 1 pixel
    subnum = 1;
    for subj = 1:5
        for subi = 1:5
            patch_sub = slide_Data_5x(strcmpi(slide_Data_5x(:,6),['sub' int2str(subnum)]),:);
            if(~isempty(patch_sub))
                out_map_5x(400*(subi-1)+1:400*subi,400*(subj-1)+1:400*subj)=out_map_5x(400*(subi-1)+1:400*subi,400*(subj-1)+1:400*subj) + patch_sub{1,3};
                count_map_5x(400*(subi-1)+1:400*subi,400*(subj-1)+1:400*subj)=count_map_5x(400*(subi-1)+1:400*subi,400*(subj-1)+1:400*subj) + 1;
                %test_mask = zeros(imgsize);
                %test_mask(400*(subi-1)+1:400*subi,400*(subj-1)+1:400*subj)=1;
                %imtool3D(bw_img,[],[],[],[],test_mask);
            end
            xsub_list = slide_Data_10x(strcmpi(slide_Data_10x(:,6),['sub' int2str(subnum)]),:);
            sub_20_pull = slide_Data_20x(strcmpi(slide_Data_20x(:,6),['sub' int2str(subnum)]),:);
            xsubnum=1;
            for xsubj=1:2
                for xsubi=1:2
                    patch_xsub = xsub_list(strcmpi(xsub_list(:,7),['xsub' int2str(xsubnum)]),:);
                    if(~isempty(patch_xsub))
                          out_map_10x(400*(subi-1)+200*(xsubi-1)+1:400*(subi-1)+200*xsubi,400*(subj-1)+200*(xsubj-1)+1:400*(subj-1)+200*xsubj)=  out_map_10x(400*(subi-1)+200*(xsubi-1)+1:400*(subi-1)+200*xsubi,400*(subj-1)+200*(xsubj-1)+1:400*(subj-1)+200*xsubj) + patch_xsub{1,3};
                        count_map_10x(400*(subi-1)+200*(xsubi-1)+1:400*(subi-1)+200*xsubi,400*(subj-1)+200*(xsubj-1)+1:400*(subj-1)+200*xsubj)=count_map_10x(400*(subi-1)+200*(xsubi-1)+1:400*(subi-1)+200*xsubi,400*(subj-1)+200*(xsubj-1)+1:400*(subj-1)+200*xsubj) + 1;
                    end
                    xxsub_list = sub_20_pull(strcmpi(sub_20_pull(:,7),['xsub' int2str(xsubnum)]),:);
                    xxsubnum = 1;
                    for xxsubj = 1:2
                        for xxsubi=1:2
                            patch_xxsub = xxsub_list(strcmpi(xxsub_list(:,8),['xxsub' int2str(xxsubnum)]),:);
                            if(~isempty(patch_xxsub))
                                  out_map_20x(400*(subi-1)+200*(xsubi-1)+100*(xxsubi-1)+1:400*(subi-1)+200*(xsubi-1)+100*xxsubi,400*(subj-1)+200*(xsubj-1)+100*(xxsubj-1)+1:400*(subj-1)+200*(xsubj-1)+100*xxsubj)=  out_map_20x(400*(subi-1)+200*(xsubi-1)+100*(xxsubi-1)+1:400*(subi-1)+200*(xsubi-1)+100*xxsubi,400*(subj-1)+200*(xsubj-1)+100*(xxsubj-1)+1:400*(subj-1)+200*(xsubj-1)+100*xxsubj) + patch_xxsub{1,3};
                                count_map_20x(400*(subi-1)+200*(xsubi-1)+100*(xxsubi-1)+1:400*(subi-1)+200*(xsubi-1)+100*xxsubi,400*(subj-1)+200*(xsubj-1)+100*(xxsubj-1)+1:400*(subj-1)+200*(xsubj-1)+100*xxsubj)=count_map_20x(400*(subi-1)+200*(xsubi-1)+100*(xxsubi-1)+1:400*(subi-1)+200*(xsubi-1)+100*xxsubi,400*(subj-1)+200*(xsubj-1)+100*(xxsubj-1)+1:400*(subj-1)+200*(xsubj-1)+100*xxsubj) + 1;
                            end
                            xxsubnum=xxsubnum+1;
                        end
                    end
                    xsubnum = xsubnum+1;
                end
            end
            subnum = subnum + 1;
        end
    end
    
    inds5 = find(count_map_5x > 0);
    out_map_5x(inds5) = out_map_5x(inds5)./count_map_5x(inds5);
    out_map_5x = out_map_5x.*tissue;
    inds10 = find(count_map_10x > 0);
    out_map_10x(inds10) = out_map_10x(inds10)./count_map_10x(inds10);
    out_map_10x = out_map_10x.*tissue;
    inds20 = find(count_map_20x > 0);
    out_map_20x(inds20) = out_map_20x(inds20)./count_map_20x(inds20);
    out_map_20x = out_map_20x.*tissue;
    
    
    mapped_image = (out_map_5x + out_map_10x + out_map_20x)./3;

    % Create indexed image, explicitly using 256 colors
    imInd=gray2ind(mapped_image,255);
    % Convert indexed image to RGB using 256-colors jet map
    jetRGB=ind2rgb(imInd,jet(255));

    %write to jpeg
    imwrite(out_map_5x,[mapLoc filesep tma_id  '_5x_prob_' tma_outcome '.jpeg']);
    imwrite(out_map_10x,[mapLoc filesep tma_id '_10x_prob_' tma_outcome '.jpeg']);
    imwrite(out_map_20x,[mapLoc filesep tma_id '_20x_prob_' tma_outcome '.jpeg']);
    imwrite(mapped_image,[mapLoc filesep tma_id '_sum_prob_' tma_outcome '.jpeg']);
    
    %keep as mat file
    save([mapLoc2 filesep tma_id  '_5x_prob_' tma_outcome '.mat'],'out_map_5x');
    save([mapLoc2 filesep tma_id '_10x_prob_' tma_outcome '.mat'],'out_map_10x');
    save([mapLoc2 filesep tma_id '_20x_prob_' tma_outcome '.mat'],'out_map_20x');
    save([mapLoc2 filesep tma_id '_sum_prob_' tma_outcome '.mat'],'mapped_image');
%     
    
    imwrite(jetRGB,[mapLoc filesep tma_id '_color_prob_' tma_outcome '.jpeg']);
    imwrite(tma_img,[mapLoc filesep tma_id '_HE_' tma_outcome '.png'])
    %copyfile([mainDir filesep tcga_id filesep 'voi' filesep strrep(svs_file,'.svs','.png')],[mapLoc filesep slides_unique{casei} '_HE_' patient_outcome '.png'])
    
end
    
