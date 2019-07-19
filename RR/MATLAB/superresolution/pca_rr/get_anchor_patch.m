function x_patch = get_anchor_patch(LR_LF,n_patch_strt,n_patch_end, m_patch_strt,m_patch_end, anchor_sai_idx)
if m_patch_end <= size(LR_LF,2) && n_patch_end <= size(LR_LF,1)
    x_patch = LR_LF(n_patch_strt:n_patch_end, m_patch_strt:m_patch_end,anchor_sai_idx);
elseif m_patch_end <= size(LR_LF,2) && n_patch_end > size(LR_LF,1)
    x_patch = padarray(LR_LF(n_patch_strt:size(LR_LF,1), m_patch_strt:m_patch_end,anchor_sai_idx),[n_patch_end-size(LR_LF,1),0],'replicate','post');
elseif m_patch_end > size(LR_LF,2) && n_patch_end > size(LR_LF,1)
    x_patch = padarray(LR_LF(n_patch_strt:size(LR_LF,1), m_patch_strt:size(LR_LF,2),anchor_sai_idx),[n_patch_end-size(LR_LF,1),m_patch_end-size(LR_LF,2)],'replicate','post');
elseif m_patch_end > size(LR_LF,2) && n_patch_end <= size(LR_LF,1)
    x_patch = padarray(LR_LF(n_patch_strt:n_patch_end, m_patch_strt:size(LR_LF,2),anchor_sai_idx),[0, m_patch_end-size(LR_LF,2)],'replicate','post');
end