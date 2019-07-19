function window_i = get_search_window_sub_ap_img(LR_LF,search_window,idx, patch_size)
window_i = LR_LF(search_window.n_strt:search_window.n_end, search_window.m_strt:search_window.m_end,idx);

if size(window_i,1) < patch_size && size(window_i,2) < patch_size
    window_i = padarray(window_i,[patch_size - size(window_i,1), patch_size - size(window_i,2)],'replicate','post');
elseif size(window_i,1) < patch_size
    window_i = padarray(window_i,[patch_size - size(window_i,1),0],'replicate','post');
elseif size(window_i,2) < patch_size
    window_i = padarray(window_i,[0, patch_size - size(window_i,2)],'replicate','post');
end