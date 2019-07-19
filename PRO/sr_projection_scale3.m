function img_output = sr_projection_scale3(lf_gray_ds, disparity, ds_flag)
% This function is used to super resolute the central view of a light field
% dataset with the downsampled grayscale LF data and the central view's
% unit disparity.

%% generate disparity and downsample the data
scale = 3;
[U,V,X,Y] = size(lf_gray_ds);%get the parameters of donwsampled ligth field
Cv = ceil(U/2);
disparity = generate_disparity_map_estimate(disparity,U,V);

%% super resolution by disparity map
threshold_lumin = 0.197;

disp(strcat('threshold of luminance is ',num2str(threshold_lumin),'-------'));
disp(' ');
Mesh_X = zeros([U*V*X*Y,1]);
Mesh_Y = zeros([U*V*X*Y,1]);
Mesh_V = zeros([U*V*X*Y,1]);%allocate space for Mesh

round_para =3;

srim = zeros([X*scale,Y*scale,1],'single');%the buffer to store the super-resolved image

ori_ds = squeeze(lf_gray_ds(Cv,Cv,:,:));
inter_img = imresize(ori_ds,scale);

index = 1;
for u = 1:U
    for v = 1:V
        disp(strcat('View -u:',num2str(u),' -v:',num2str(v),' is meshing into griddata.....'));
        sub_img = squeeze(lf_gray_ds(u,v,:,:));%get sub_image
        sub_disparity_x = squeeze(disparity(u,v,:,:,1));
        sub_disparity_y = squeeze(disparity(u,v,:,:,2));%get disparity

        for x = 1:X
            for y = 1:Y
                shift_x = sub_disparity_x(x,y);
                shift_y = sub_disparity_y(x,y);
                inter_index_x = x-shift_x;
                inter_index_y = y+shift_y;%x-shift,y+shift,used for HCI dataset
                

                md_x = round(inter_index_x*round_para)-1;
                md_y = round(inter_index_y*round_para)-1;
                              
                if (md_x>0) && (md_x <= X*round_para) && (md_y>0) && (md_y <= Y*round_para)
%                     if(abs(double(sub_img(x,y))-double(inter_img(md_x,md_y)))<double(inter_img(md_x,md_y))*threshold_lumin)
                    if(abs(double(sub_img(x,y))-double(inter_img(md_x,md_y)))<double(inter_img(md_x,md_y))*threshold_lumin)
                        Mesh_X(index,1) = inter_index_x;
                        Mesh_Y(index,1) = inter_index_y;
                        Mesh_V(index,1) = double(sub_img(x,y));
                        index = index+1;
                    end
                else
                    Mesh_X(index,1) = inter_index_x;
                    Mesh_Y(index,1) = inter_index_y;
                    Mesh_V(index,1) = double(sub_img(x,y));
                    index = index+1;
                end
                %----------------------------------------------
                
            end
        end
    end
end


%% delete the zero points at the tail of Mesh vectors(FOR THE ONCE DELETE ALGO)
Mesh_X = Mesh_X(1:index-1,1);
Mesh_Y = Mesh_Y(1:index-1,1);
Mesh_V = Mesh_V(1:index-1,1);
del_percent = 1.0 - (index-1)/(X*Y*U*V);
disp(strcat('The proportion of deleting points is:',num2str(del_percent)));

%% interpolation and sampling
[Mesh_Xi,Mesh_Yi] = meshgrid(2.0/3.0:1.0/3.0:Y+1.0/3.0,2.0/3.0:1.0/3.0:X+1.0/3.0);
disp('Mesh done!');
disp(' ');
disp('Start Interpolation');
% tic
Surf = scatteredInterpolant(Mesh_Y,Mesh_X,Mesh_V);

num_temp = Surf(Mesh_Xi,Mesh_Yi);
srim(:,:,1) = round(num_temp);
srim = uint8(srim);

%% iteration to get the best one 
imgH0 = inter_img;
imgH1 = srim;
imgL0 = ori_ds;
% iteration 1
imgH2 = iteration_once(imgH1,imgL0,scale,ds_flag);
psnr2 = psnr(imgH0,imgH2);
% iteration 2
imgH3 = iteration_once(imgH2,imgL0,scale,ds_flag);
psnr3 = psnr(imgH0,imgH3);
iter_num = 2;
% more iterations
while abs(psnr3-psnr2)>0
    imgH3 = iteration_once(imgH3,imgL0,scale,ds_flag);
    iter_num = iter_num + 1;
    psnr2 = psnr3;
    psnr3 = psnr(imgH0,imgH3);
end

img_output = imgH3;



