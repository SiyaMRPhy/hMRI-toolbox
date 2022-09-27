function ss_eTPM_masked()
% genetrate TPM maksed (no eyes), c456 combined
% this TPM required to estimate the masked PD values  
%--------------------------------------------------------------------------
% SIYA SHERIF
% 15JUL2022
% CRC - ULiege, Liege

close all; clear all; clc;

% path to the eTPM  (i am runing it in tets folder, NOT touching the
% original TPM folder, will copy the scipt and output to the hMRI TPM folder);
eTPM = fullfile(pwd,'eTPM.nii') ;

% tpm_out_name = spm_file(eTPM_tmp,'suffix','_mask');
ss_eTPM_c456 = spm_file(eTPM,'prefix','ss_','suffix','_c456');

if ~exist(ss_eTPM_c456)
    
    % get eTPM
    Y_eTPM =  spm_vol(eTPM);

    V_c1 =  double(spm_read_vols(Y_eTPM(1)));
    V_c2 =  double(spm_read_vols(Y_eTPM(2)));
    V_c3 =  double(spm_read_vols(Y_eTPM(3)));
    
    % ICV mask
    maskICV = fullfile(spm_file(eTPM,'path'),'maskICV_eTPM.nii');
    Y_maskICV=  spm_vol(maskICV);
    V_maskICV =  double(spm_read_vols(Y_maskICV));
    
    ss_eTPM_c1 = spm_file(eTPM,'prefix','ss_','suffix','_c1');
    if ~exist(ss_eTPM_c1)

        V_c1_tmp = (V_c1.* V_maskICV) + (~V_maskICV.*1e-6);
        clear OP
        OP          = Y_eTPM(1);
        OP.fname    = ss_eTPM_c1;%
        OP.descrip  = sprintf('modified_tpm_c1');
        spm_write_vol(OP,V_c1_tmp);
        
    end
    
    ss_eTPM_c2 = spm_file(eTPM,'prefix','ss_','suffix','_c2');
    if ~exist(ss_eTPM_c2)

        V_c2_tmp = (V_c2.*V_maskICV) + (~V_maskICV.*1e-6);
        clear OP
        OP          = Y_eTPM(1);
        OP.fname    = ss_eTPM_c2;%
        OP.descrip  = sprintf('modified_tpm_c2');
        spm_write_vol(OP,V_c2_tmp);
        
    end
    
    ss_eTPM_c3 = spm_file(eTPM,'prefix','ss_','suffix','_c3');
    if ~exist(ss_eTPM_c3)
        
        V_c3_tmp = (V_c3.*V_maskICV) + (~V_maskICV.*1e-6);
        % mask and zeros to 10e^-6
        clear OP
        OP          = Y_eTPM(1);
        OP.fname    = ss_eTPM_c3;%
        OP.descrip  = sprintf('modified_tpm_c3');
        spm_write_vol(OP,V_c3_tmp);
        
    end

    ss_eTPM_c456 = spm_file(eTPM,'prefix','ss_','suffix','_c456');
    
    if ~exist(ss_eTPM_c456)
        
        
        V_c456_tmp = 1-(V_c1_tmp + V_c2_tmp + V_c3_tmp); % sum of all probability should add to 1 
        clear OP
        OP          = Y_eTPM(1);
        OP.fname    = ss_eTPM_c456;%
        OP.descrip  = sprintf('modified_tpm_c456');
        spm_write_vol(OP,V_c456_tmp);
        
    end

end