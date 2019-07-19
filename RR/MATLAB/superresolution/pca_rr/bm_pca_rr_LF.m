function SR_LF = bm_pca_rr_LF(LR_LF,param)
% This script is used to run the block-matck pca ridge-regression algorithm. 
% This function receives the low resolution lightfield LR_LF and super-resolved
% it to a resolution defined by sz. This method divides the lightfield into 4D patch
% volumes, were the patches from each sub-aparture image are selected using block
% matching since this will align the lightfield and thus reduce the disparity. 
% Reducing the disparity is important since it reduces the variance across the 
% patch-volume and thus makes it easier for the super-resolution algorithm to find
% a good upscaling projection matrix.
%
% Input: LR_LF: This corresponds to the low-resolution lightfield
%        param: structure containing different parameters, dictionaries and pca 
%               basis functions which were pre-computed
%        sz:    Size of the target high-resolution lightfield 
%
% Output: SR_LF:The super-resolved lightfield
% 
% Reuben Farrugia
% version 1: 23/8/2016
% version 2: 26/1/2018
%

% Compute bicubic interpolation for each sub-aparture image separately
%LR_LF = bicubic_LF_upscale(LR_LF,sz);

% Initialize the super-resolved lightfield
SR_LF = -ones(size(LR_LF)); % -1 indicates that the pixel is empty and needs to
                            % be approximated. All pixels in SR_LF are set to -1

%--------------------------------------------------------------------------
% Configuration
%--------------------------------------------------------------------------
patch_size = param.patch_size;  % patchsize in the xy dimension
overlap    = param.overlap;     % overlap in the xy dimension
window_size = param.window_size;% Window size

% Extract information from the pca basis
El = param.pca_basis.El; Ml = param.pca_basis.Ml; 
Eh = param.pca_basis.Eh; Mh = param.pca_basis.Mh; 

% Get the projection matrix for ridge regression
phi = param.phi;
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Load information needed to compute LF SR patch-based
%-------------------------------------------------------------------------------

% Derive the dimensions of the image
H = size(LR_LF,1); W = size(LR_LF,2);

% Derive the number of horizontal patches
Nh = ceil(W/(patch_size-overlap));

% Derive the number of vertical patches
Nv = ceil(H/(patch_size-overlap));

% Derive the index of sub-aparture images having holes
idx_holes = 1:size(SR_LF,3); 

while ~isempty(idx_holes)
    if size(idx_holes,2) == size(LR_LF,3)
        % Set the anchor_sai_idx to the center sub-aparture image
        anchor_sai_idx = median(1:size(LR_LF,3));
        % Find the index of the ahchor view
        idx = find(idx_holes == anchor_sai_idx);
    else
        % Determine the number of holes for each sub-aparture image
        Nholes = reshape(sum(sum(SR_LF == -1,1),2),[size(SR_LF,3),1]);
    
        % Consider only those marked in idx_holes
        Nholes = Nholes(idx_holes);
    
        % Derive the sub-aparture image with the smallest number of holes
        [~,idx] = min(Nholes);
    
        % Derive the sub-aparture index to be used as anchor
        anchor_sai_idx = idx_holes(idx);
    end
    
    % Initialize the n patch index
    n_patch_strt = 1;
    
    for n = 1:Nv
        % Initialize the m patch start to 1
        m_patch_strt = 1;
        for m = 1:Nh
            % Derive the m patch index
            m_patch_end = m_patch_strt + patch_size - 1;
            % Derive the n patch index
            n_patch_end = n_patch_strt + patch_size - 1;
      
            % There are some pixels that need to be restored
            if m_patch_end <= size(LR_LF,2) && n_patch_end <= size(LR_LF,1)
                x_patch = LR_LF(n_patch_strt:n_patch_end, m_patch_strt:m_patch_end,anchor_sai_idx);
                mask    = SR_LF(n_patch_strt:n_patch_end, m_patch_strt:m_patch_end,anchor_sai_idx) == -1;
            elseif m_patch_end <= size(LR_LF,2) && n_patch_end > size(LR_LF,1)
                x_patch = padarray(LR_LF(n_patch_strt:size(LR_LF,1), m_patch_strt:m_patch_end,anchor_sai_idx),[n_patch_end-size(LR_LF,1),0],'replicate','post');
                mask    = padarray(SR_LF(n_patch_strt:size(SR_LF,1), m_patch_strt:m_patch_end,anchor_sai_idx),[n_patch_end-size(SR_LF,1),0],'replicate','post') == -1;
            elseif m_patch_end > size(LR_LF,2) && n_patch_end > size(LR_LF,1)
                x_patch = padarray(LR_LF(n_patch_strt:size(LR_LF,1), m_patch_strt:size(LR_LF,2),anchor_sai_idx),[n_patch_end-size(LR_LF,1),m_patch_end-size(LR_LF,2)],'replicate','post');
                mask    = padarray(SR_LF(n_patch_strt:size(SR_LF,1), m_patch_strt:size(SR_LF,2),anchor_sai_idx),[n_patch_end-size(SR_LF,1),m_patch_end-size(SR_LF,2)],'replicate','post') == -1;
            elseif m_patch_end > size(LR_LF,2) && n_patch_end <= size(LR_LF,1)
                x_patch = padarray(LR_LF(n_patch_strt:n_patch_end, m_patch_strt:size(LR_LF,2),anchor_sai_idx),[0, m_patch_end-size(LR_LF,2)],'replicate','post');
                mask    = padarray(SR_LF(n_patch_strt:n_patch_end, m_patch_strt:size(SR_LF,2),anchor_sai_idx),[0, m_patch_end-size(SR_LF,2)],'replicate','post') == -1;
            end
            
            if sum(sum(mask)) > 0 % This patch has to be restored            
                               
                % Derive the search window
                search_window = get_search_window(n_patch_strt, n_patch_end, m_patch_strt, m_patch_end,size(LR_LF),window_size);
    
                % Derive the low-resolution aligned patch volume
                [LR_patch_volume,dv] = get_LR_aligned_patch_volume(LR_LF,search_window,x_patch,patch_size,anchor_sai_idx);
        
                %--------------------------------------------------------------------------
                % Super-resolve the volume patch
                %--------------------------------------------------------------------------
                % Extract the low-quality patch
                xp = double(reshape(LR_patch_volume,[patch_size*patch_size*size(LR_LF,3),1]));
        
                % Center the low-resolution patch
                xp_c = xp - Ml;
        
                % Project the low-resolution sample on the subjspace
                xp_ss = El' * xp_c;
    
                % Approximate the location on high-resolution subspace
                yp_ss_hat = phi * xp_ss;
    
                % Project back to get the high-resolution patch
                yp_c_hat = Eh * yp_ss_hat;
            
                % Add the mean 
                yp_hat = round(yp_c_hat + Mh);

                % Reshape the high-quality volume patch
                yp = reshape(yp_hat,[patch_size,patch_size,size(LR_LF,3)]);
       
                %--------------------------------------------------------------------------
                % Reconstruct the high-resolution lightfield
                %--------------------------------------------------------------------------
                for i = 1:size(SR_LF,3)
                    if i == anchor_sai_idx
                        idx_m = m_patch_strt:min([m_patch_end,size(SR_LF,2)]);
                        idx_n = n_patch_strt:min([n_patch_end,size(SR_LF,1)]);
                        xp_i = SR_LF(idx_n,idx_m,i);
                        mask_i = ~(xp_i == -1);
                        % Stitch the recovered patch together
                        yp_i = stitching(xp_i,mask_i,yp(:,:,i));
                        % Update the super-resolved lightfield
                        SR_LF(idx_n, idx_m,i) = yp_i;
                    else
                        idx_m = search_window{1}.m_strt + dv(i,2) - 1:min([search_window{1}.m_strt+patch_size-1+dv(i,2)-1,W]);
                        idx_n = search_window{1}.n_strt + dv(i,1) - 1:min([search_window{1}.n_strt+patch_size-1+dv(i,1)-1,H]);
                        xp_i = SR_LF(idx_n, idx_m,i);
                        mask_i = ~(xp_i == -1);
                        % Stitch the recovered patch together
                        yp_i = stitching(xp_i,mask_i,yp(:,:,i));
                        % Update the super-resolved lightfield
                        SR_LF(idx_n, idx_m,i) = yp_i;
                    end
                end
            end    
            % Update the start values
            m_patch_strt = m_patch_end - overlap + 1;
       
        end
        % Update the n patch strt
        n_patch_strt = n_patch_end - overlap + 1;    
    end
    % Remove this index from idx_holes
    idx_holes(idx) = [];
end

% Make sure the super-resolved lightfield is integer
SR_LF = uint8(SR_LF);
