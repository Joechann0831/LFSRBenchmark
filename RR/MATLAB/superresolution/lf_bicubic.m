function SR_LF = lf_bicubic(LR_LF,sz,mf)

% Initialize the super-resolved light field
SR_LF = zeros(sz,'uint8');

for u = 1:size(LR_LF,4)
    for v = 1:size(LR_LF,5)
        % Extract the current sub-aperture image
        Ilr = LR_LF(:,:,:,u,v);
        
        Isr = imresize(Ilr,mf,'bicubic');
        
        Isr = Isr(1:sz(1),1:sz(2),:);
        
        % Put the restored sub-aperture image in the restored light field
        SR_LF(:,:,:,u,v) = Isr;
    end
end