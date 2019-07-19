function depth_output_n = DEPTH_NORMALIZE(depth_output)
depth_output_min = min(min(depth_output))                                 ;
depth_output_max = max(max(depth_output))                                 ;
depth_output_n    = (depth_output-depth_output_min)/...
                                       (depth_output_max-depth_output_min);
end