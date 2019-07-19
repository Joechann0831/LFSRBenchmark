function LF_pvs = lf2patchvol(LF, patch_size, overlap)
% This function divides a lightfield into overlapping patch volues of size
% patch_size x patch_size x #sub-aparture images. This function recenives a
% lightfields and returns a list of overlapping patch volumes
%
% Input: LF - lightfield tensor
%        patch_size - size of the patch in the xy dimension
%        overlap - overlap of the patch volume
%
% Output: LF_pvs - list of patch volumes extracted from the lightifeld
%
% Reuben Farrugia
% Date: 30/07/2016

% Derive the dimensions of the image
H = size(LF,1); W = size(LF,2);

% Derive the number of horizontal patches
Nh = ceil(W/(patch_size-overlap));

% Derive the number of vertical patches
Nv = ceil(H/(patch_size-overlap));

% Derive the number of patches
N = Nh * Nv;

% Initialize x to contain N patches
LF_pvs = cell(N,1);

% Initialize the n patch index
n_patch_strt = 1;
% Initialize the patch index
k = 1;
for n = 1:Nv
    % Initialize the m patch start to 1
    m_patch_strt = 1;
    for m = 1:Nh
        % Derive the m patch index
        m_patch_end = m_patch_strt + patch_size - 1;
        % Derive the n patch index
        n_patch_end = n_patch_strt + patch_size - 1;
        
        % Extract the kth patch
        if m_patch_end <= W && n_patch_end <= H
            x = LF(n_patch_strt:n_patch_end, m_patch_strt:m_patch_end,:);
        elseif m_patch_end <= W && n_patch_end > H
            x = padarray(LF(n_patch_strt:H, m_patch_strt:m_patch_end,:),[n_patch_end-H,0,0],'replicate','post');
            %[LF(n_patch_strt:H, m_patch_strt:m_patch_end,:); zeros(n_patch_end-H,patch_size,size(LF,3))];
        elseif m_patch_end > W && n_patch_end > H
            x = padarray(LF(n_patch_strt:H, m_patch_strt:W,:),[n_patch_end-H,m_patch_end-W,0],'replicate','post');
            %x = [[LF(n_patch_strt:H, m_patch_strt:W,:), zeros(H-n_patch_strt+1,m_patch_end-W,size(LF,3))]; zeros(n_patch_end-H,patch_size,size(LF,3))];
        elseif m_patch_end > W && n_patch_end <= H
            x = padarray(LF(n_patch_strt:n_patch_end, m_patch_strt:W,:),[0, m_patch_end-W,0],'replicate','post');
            %x = [LF(n_patch_strt:n_patch_end, m_patch_strt:W,:), zeros(patch_size,m_patch_end-W,size(LF,3))];
        end
        % Store the vectorized representation of the volume patch
        LF_pvs{k,1} = reshape(x,[patch_size*patch_size*size(LF,3)*size(LF,4),1]);
        
        % Update the start values
        m_patch_strt = m_patch_end - overlap + 1;
        
        % Increment the patch counter
        k = k + 1;
    end
    % Update the n patch strt
    n_patch_strt = n_patch_end - overlap + 1;
end
