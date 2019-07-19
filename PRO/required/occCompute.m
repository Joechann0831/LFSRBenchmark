function occlusion = occCompute(d_var, d_foc, im_edge, dir, v_max, f_max, occ_thre)

% correspondence cue
occ_v = max(d_var(:,:,1)./d_var(:,:,2), d_var(:,:,2)./d_var(:,:,1)) .* (dir>-100); 
occ_v(~isfinite(occ_v)) = 0; 
occ_v(occ_v > v_max) = v_max;
occ_v = (occ_v-mean(occ_v(:))) / std(occ_v(:));
% refocus cue
occ_f = sum((abs(d_foc(:,:,:,1)-d_foc(:,:,:,2))).^2, 3) .* (dir>-100);
occ_f(~isfinite(occ_f)) = 0; 
occ_f(occ_f > f_max) = f_max;
occ_f = (occ_f-mean(occ_f(:))) / std(occ_f(:));
% overall response
occ = bwareaopen(occ_v.*occ_f > 1, 2);
occlusion = occEdge(occ, im_edge, occ_thre);