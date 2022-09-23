function ss_maskICV_eTPM()
% generate a ICV mask for eTPM 
% MASK is cerated by thresholding the eTPM
% this required to estimate the masked PD values  
% uses fucntion fill_it, connect_it , open_it from the FieldMap toolbox
%--------------------------------------------------------------------------
% SIYA SHERIF
% 15JUL2022
% CRC - ULiege, Liege
% should we use this for the estimataing the c1,c2,c3 values ? ask chris or SPM group

close all; clear all; clc;

% path to the eTPM  (i am runing it in tets folder, NOT touching the
% original TPM folder, will copy the scipt and output to the hMRI TPM folder);
eTPM = fullfile(pwd,'eTPM.nii') ;

% define output name
maskICV =  spm_file(eTPM,'prefix','maskICV_');

if ~exist(maskICV)
    
    % read the eTPM vol
    Y_eTPM =  spm_vol(eTPM);
    V_c1 =  double(spm_read_vols(Y_eTPM(1)));
    V_c2 =  double(spm_read_vols(Y_eTPM(2)));
    V_c3 =  double(spm_read_vols(Y_eTPM(3)));
    V_c4 =  double(spm_read_vols(Y_eTPM(4)));
    V_c5 =  double(spm_read_vols(Y_eTPM(5)));
    V_c6 =  double(spm_read_vols(Y_eTPM(6)));
    
    % create a mask by thresholding  c1, c2, c3
    c_thresh = 0.2;
    
    % threshold each tissue class separately and combine
    TPM_mask = ((V_c1 > c_thresh)  + (V_c2 > c_thresh) + (V_c3 > c_thresh)) > 0 ;
    
    % all other issue class combined 
    notICV =  ((V_c4 + V_c5 + V_c6) > 0.6);
    
    % multiply with the ~(non ICV)
    TPM_mask = TPM_mask .* (~(notICV));
    
    % this part is from the SPM Fieldmap 
    nerode  = 2;
    ndilate = 2;
    thresh  = 0.5;  
    fwhm    = 5;   
    
    TPM_mask = open_it(TPM_mask,nerode,ndilate); % Do opening to get rid of scalp

    % Calculate kernel in voxels:
    vxs = sqrt(sum(Y_eTPM(1).mat(1:3,1:3).^2));
    fwhm = repmat(fwhm,1,3)./vxs;
    TPM_mask = fill_it(TPM_mask,fwhm,thresh); % Do fill to fill holes
    
    % removing the 5 lines in the z plane (simillar to the TPM mask in SPM )
    % need to ask the logic behind removing 5 line in this plane to chris or SPM group 
    TPM_mask(:,:,1:5) = 0; 

    OP          = Y_eTPM(1);
    OP.fname    = maskICV;%
    OP.descrip  = sprintf('Mask:erode=%d,dilate=%d,fwhm=%d,thresh=%1.1f',nerode,ndilate,fwhm,thresh);
    spm_write_vol(OP,TPM_mask);

end



function ovol=open_it(vol,ne,nd)
    % Do a morphological opening. This consists of an erosion, followed by 
    % finding the largest connected component, followed by a dilation.

    % Do an erosion then a connected components then a dilation 
    % to get rid of stuff outside brain.
    for i=1:ne
       nvol=spm_erode(double(vol));
       vol=nvol;
    end
    nvol=connect_it(vol);
    vol=nvol;
    for i=1:nd
       nvol=spm_dilate(double(vol));
       vol=nvol;
    end

    ovol=nvol;
end




function ovol=fill_it(vol,k,thresh)
    % Do morpholigical fill. This consists of finding the largest connected 
    % component and assuming that is outside of the head. All the other 
    % components are set to 1 (in the mask). The result is then smoothed by k
    % and thresholded by thresh.
    ovol=vol;

    % Need to find connected components of negative volume
    vol=~vol;
    [vol,NUM]=spm_bwlabel(double(vol),26); 

    % Now get biggest component and assume this is outside head..
    pnc=0;
    maxnc=1;
    for i=1:NUM
       nc=size(find(vol==i),1);
       if nc>pnc
          maxnc=i;
          pnc=nc;
       end
    end

    % We know maxnc is largest cluster outside brain, so lets make all the
    % others = 1.
    for i=1:NUM
        if i~=maxnc
           ovol(vol==i)=1;
        end
    end

    spm_smooth(ovol,ovol,k);
    ovol=ovol>thresh;
end

function ovol=connect_it(vol)
    % Find connected components and return the largest one.

    [vol,NUM]=spm_bwlabel(double(vol),26); 

    % Get biggest component
    pnc=0;
    maxnc=1;
    for i=1:NUM
       nc=size(find(vol==i),1);
       if nc>pnc
          maxnc=i;
          pnc=nc;
       end
    end
    ovol=(vol==maxnc);
end

end