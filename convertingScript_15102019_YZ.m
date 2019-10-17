clear
load('E:\Yue_folder\OneDrive - bwstaff\AA LAB\Miniprojects-master\python files\example01.mat')
[real_sphere,real_status,rearranged] = match_led(vertices,stimulus');
folder_name = 'E:\Yue_folder\OneDrive - bwstaff\AA LAB\Miniprojects-master\Output files\';
sti_pattern.folder_name     = folder_name;
sti_pattern.Pattername      = {['Left_',datestr(datetime,30),'_Pattern'];['Right_',datestr(datetime,30),'_Pattern']};
sti_pattern.x_num           = size(stimulus,1);
sti_pattern.y_num           = 1;
sti_pattern.num_panels      = 120;
sti_pattern.gs_val          = 2;
sti_pattern.row_compression = 0;
sti_pattern.Panel_map       = flipud(reshape(1:120,8,15));
sti_pattern.CFname          = {'CF1.img','CF2.img'};
for lr = 1:2
    temp_pat                = sti_pattern;
    temp_pat.Pats           = double(rearranged(:,1+(end/2)*(lr-1):(end/2)*lr,:));
    temp_pat.Pats(isnan(temp_pat.Pats)) = 0;
    temp_pat.BitMapIndex    = process_panel_map(temp_pat);
    temp_pat.data           = Make_pattern_vector(temp_pat,1);
    
    
    %%% Part Below for making CF images
    current_frame_size      = temp_pat.num_panels*temp_pat.gs_val*8;
    current_num_frames      = temp_pat.x_num*temp_pat.y_num;
    
    block_size              = 512;
    block_offset            = 128;
    block_indexer           = 1;
    blocks_per_frame        = ceil(current_frame_size/block_size);
    block_start_address     = block_offset + block_indexer;
    
    Header_block            = zeros(1, block_size);
    Header_block(1:6)       = [dec2char(temp_pat.x_num,2), dec2char(temp_pat.y_num,2), temp_pat.num_panels, temp_pat.gs_val];
    Header_block(7:10)      = dec2char(block_start_address,4);
    
    % now write all of the frame info
    for i = 1:current_num_frames
        cf_start_address = (block_indexer - 1)*block_size + 1;
        cf_end_address = cf_start_address + current_frame_size - 1;
        % always forced to start frame at a block boundary
        pat_start_address = (i - 1)*current_frame_size + 1;
        pat_end_address = pat_start_address + current_frame_size - 1;
        Pattern_Data(cf_start_address:cf_end_address) = temp_pat.data(pat_start_address:pat_end_address);
        block_indexer = block_indexer + blocks_per_frame;
    end
    
    
    Save CF img
    fid = fopen([sti_pattern.folder_name,sti_pattern.CFname{lr}], 'w');
    fwrite(fid, Header_block,'uchar');
    fclose(fid);
    fid = fopen([sti_pattern.folder_name,sti_pattern.CFname{lr}], 'a');
    fwrite(fid, Pattern_Data,'uchar');
    fclose(fid);
    % Save matlab Pattern
    save([sti_pattern.folder_name,sti_pattern.Pattername{lr},'.mat'],'temp_pat');
end
%% Demo video of the stimulus
clf
[vx,vy,vz] = sphere(30);
surf(vx.*100,vy.*100,vz.*100,'facecolor',ones(1,3)./2,'edgecolor','none')
hold on
demos = scatter3(real_sphere(:,1),real_sphere(:,2),real_sphere(:,3),10,real_status(:,1).*ones(1,3),'filled');
axis image off
colormap gray
%
for i = 1:1000
    demos.CData = real_status(:,i).*ones(1,3);
    pause(.01)
end