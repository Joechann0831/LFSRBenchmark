function search_window = get_search_window(n_patch_strt, n_patch_end, m_patch_strt, m_patch_end,sz,delta)
% This function returns the search window for each sup-aparture image where
% the window-size is dependent on the disparity
search_window = cell(size(delta,1),1);
for i = 1:size(delta,1)
    search_window{i}.n_strt = max([n_patch_strt-delta, 1]);
    search_window{i}.n_end  = min([n_patch_end+delta,sz(1)]);
    search_window{i}.m_strt = max([m_patch_strt-delta, 1]);
    search_window{i}.m_end  = min([m_patch_end+delta,sz(2)]);
end