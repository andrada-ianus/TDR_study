

% path to current directory
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);


 


%% Code to calculate TDR

% labels for the data from the two segments of the spinal cord  
% G600 refers to the data acquired with the limited gradient strength,
% while Gmax refers to the data acquired with maximum gradient strength 
% see Warner et al for the details of the acquisition


ids = {'Aeon_SC_New_2_Gmax','Aeon_SC_New_2_G600','Aeon_SC_New_1_Gmax','Aeon_SC_New_1_G600'};


dims = size(image_den);


for i = 1 
    
     mouse_id = ids{i};
     data_folder = [filepath '/data/' mouse_id];
     
     % loading the denoised data
     load([data_folder filesep mouse_id '_Den']);
   

     image_den_norm = zeros(size(image_den));
        
       %  voxelwise_sorted_indices = zeros([dims(1:3) 60]); % 10 b0s and 60 directions;
     % gtab is a structure which contains information about the acquisition
     % protocol, such as bval, bvec, bmat, diffusion time (delta) and
     % gradient duration (smalldel)
     % gtab.expno - refers to each shell
     
     
    % each shell is normalized by its respective b0 values
        for expno = 1:max(gtab.expno)
            b0_indices = find(gtab.bval <100 & gtab.expno == expno);
            meanS0 = mean(abs(image_den(:,:,:,b0_indices)),4);

            image_den_norm(:,:,:,gtab.expno == expno) = abs(image_den(:,:,:,gtab.expno == expno))./meanS0;
        end
        
        
        % sorting the measurements based on the sum of the signal values in
        % the three shells, as explained in Warner et all. For this
        % acquisition, the order of the gradient directions is the same for
        % the different shells. This can be seen in the gtab.bvec variable.
        
        [tmp, voxelwise_sorted_indices] = sort(image_den_norm(:,:,:,11:70)+image_den_norm(:,:,:,81:140)+image_den_norm(:,:,:,151:210),4,'desc');
        
       % non optimised TDR (calculated from the shells with long diffusion time and short diffusion time) 
       TDR_orig = zeros([dims(1:3) 60]); 
       
       % optimised TDR (calculated from the shells with long gradient duration and short diffusion time) 
       TDR_opt = zeros([dims(1:3) 60]); 
       
       for it1 = 1:dims(1)
           for it2 = 1:dims(2)
               for it3 = 1:dims(3)
       
                   for it4 = 1:60

                       TDR_orig(it1,it2,it3,it4) = (mean(image_den_norm(it1,it2,it3,80+voxelwise_sorted_indices(it1,it2,it3,1:it4)),4) - ...
                           mean(image_den_norm(it1,it2,it3,10+voxelwise_sorted_indices(it1,it2,it3,1:it4)),4))./mean(image_den_norm(it1,it2,it3,80+voxelwise_sorted_indices(it1,it2,it3,1:it4)),4);
                      
                       TDR_opt(it1,it2,it3,it4) = (mean(image_den_norm(it1,it2,it3,150+voxelwise_sorted_indices(it1,it2,it3,1:it4)),4) - ...
                       mean(image_den_norm(it1,it2,it3,10+voxelwise_sorted_indices(it1,it2,it3,1:it4)),4))./mean(image_den_norm(it1,it2,it3,150+voxelwise_sorted_indices(it1,it2,it3,1:it4)),4);
                   end
               end
           end
       end

       
%       save([data_folder filesep mouse_id '_TDR'],'TDR_orig','TDR_opt');
       
   
end

display('finished calculating TDR');
%%  load data and masks for one of the datasets


ids = {'Aeon_SC_New_1_Gmax','Aeon_SC_New_1_G600','Aeon_SC_New_2_Gmax','Aeon_SC_New_2_G600'};


i = 2 ;
    
mouse_id = ids{i}; 

 
data_folder = [filepath '/data/' mouse_id];
load([data_folder filesep mouse_id '_TDR']);
load([data_folder filesep mouse_id '_Den']);
load([data_folder filesep mouse_id '_mask_WM']);
load([data_folder filesep mouse_id '_mask_GM']);

%%  Plot GM and WM ROIs
     
figure(1);
imagesc(image_den(:,:,3,1).*mask(:,:,3)); colormap(gray); axis off; axis square; hold on;
contour(mask_WM(:,:,3),'Linecolor',[0 0.4470 0.7410],'Linewidth',1); axis square; hold on;
contour(mask_GM(:,:,3),'linecolor',[0.8500 0.3250 0.0980],'Linewidth',1); axis square

%% Diffusion weighted plots

% gradient parallel to fibre directions
figure(100);
% 1st shell
imagesc(image_den(:,:,3,11).*mask(:,:,3),[0 120]); colormap(gray);  axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)

figure(101);   
% 2nd shell
imagesc(image_den(:,:,3,81).*mask(:,:,3),[0 120]);colormap(gray); axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)

figure(102);     
% 3rd shell
imagesc(image_den(:,:,3,151).*mask(:,:,3),[0 120]); colormap(gray); axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)


% gradient perpendicular to fibre directions
figure(103);
imagesc(image_den(:,:,3,70).*mask(:,:,3),[0 120]); colormap(gray); axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)

figure(104);     
imagesc(image_den(:,:,3,140).*mask(:,:,3),[0 120]); colormap(gray); axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)

figure(105);     
imagesc(image_den(:,:,3,210).*mask(:,:,3),[0 120]); colormap(gray); axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)

% b0 image
figure(106);
imagesc(image_den(:,:,3,1).*mask(:,:,3),[0 500]); colormap(gray);  axis off; axis square; hold on; colorbar;  set(gca,'FontSize',18)
    
    

     
%% TDR map

figure(2);

TDR_plot = TDR_opt(:,:,:,12);
TDR_plot(mask == 0) = -2;
cmap = parula(256);
cmap(1,:) = 0;
imagesc(TDR_plot(:,:,3),[-0.4 0.6]); axis square; axis off; colormap(cmap); colorbar;
set(gca,'FontSize',18)
     
   

%% Angular plots - WM

TDR_plot_orig = TDR_orig;

figure(3);      
TDR_plot_orig(repmat(mask_WM,[1,1,1,60])==0) = NaN;
e = errorbar(squeeze(nanmean(nanmean(nanmean(TDR_plot_orig,1),2),3)),squeeze(nanstd(nanstd(nanstd(TDR_plot_orig,[],1),[],2),[],3)),'*');
xlabel('Number of averaged directions');
ylabel('TDR')
ylim([0 0.3])
%title('Original')
set(e,'Linewidth',2,'Color',[0 0.4470 0.7410])
set(gca,'FontSize',18)


TDR_plot_opt = TDR_opt;

figure(4);     
TDR_plot_opt(repmat(mask_WM,[1,1,1,60])==0) = NaN;
e = errorbar(squeeze(nanmean(nanmean(nanmean(TDR_plot_opt,1),2),3)),squeeze(nanstd(nanstd(nanstd(TDR_plot_opt,[],1),[],2),[],3)),'*')
xlabel('Number of averaged directions');
ylabel('TDR')
ylim([0 0.3])
% title('Optimised')
set(e,'Linewidth',2,'Color',[0 0.4470 0.7410])
set(gca,'FontSize',18)
    
  
    
   
    
%% Angular plots - GM

TDR_plot_orig = TDR_orig;
figure(5);      
TDR_orig(repmat(mask_GM,[1,1,1,60])==0) = NaN;
e = errorbar(squeeze(nanmean(nanmean(nanmean(TDR_orig,1),2),3)),squeeze(nanstd(nanstd(nanstd(TDR_orig,[],1),[],2),[],3)),'*');
xlabel('Number of averaged directions');
ylabel('TDR')
%ylim([-0.4 -0.1])
%title('Original')
set(e,'Linewidth',2,'Color',[0.8500 0.3250 0.0980])
set(gca,'FontSize',18)


TDR_plot_opt = TDR_opt;
figure(6);     
TDR_opt(repmat(mask_GM,[1,1,1,60])==0) = NaN;
e = errorbar(squeeze(nanmean(nanmean(nanmean(TDR_opt,1),2),3)),squeeze(nanstd(nanstd(nanstd(TDR_opt,[],1),[],2),[],3)),'*')
xlabel('Number of averaged directions');
ylabel('TDR')
%ylim([-0.4 -0.1])
% title('Optimised')
set(e,'Linewidth',2,'Color',[0.8500 0.3250 0.0980])
set(gca,'FontSize',18)
    
    
    

%% Axon diameter

rois = {'VST','ReST','RST','FC','dCST','FG'};

ROI_ax = [4.47,2.22,3.39,3.73,1.16,1.80]; % axon diameter values from literature

ROI_values = zeros(2,6);


ids = {'Aeon_SC_New_1_G600','Aeon_SC_New_2_G600'};

%ids = {'Aeon_SC_New_1_Gmax','Aeon_SC_New_2_Gmax'};

% combine the data fromt he two segments for each gradient setting
for i = 1:2
    
     mouse_id = ids{i};
     
     data_folder = [filepath '/data/' mouse_id];
     load([data_folder filesep mouse_id '_TDR']);
     load([data_folder filesep mouse_id '_Den']);
     load([data_folder filesep mouse_id '_mask_WM_ROIs']);
     load([data_folder filesep mouse_id '_mask_GM']);

    for j = 1:6 % the number of ROIs
       
        if contains(mouse_id,'_2')
            slice = 5;
        else 
            slice = 3;
        end
         mask_plot = mask_WM_ROI(:,:,slice,j);
        mask_plot(mask_plot == 0) = NaN;
        ROI_values(i,j) = nanmean(nanmean(TDR_opt(:,:,slice,12).*mask_plot));
      
       
    end
end


c = lines(6);

scatter(ROI_ax,ROI_values(1,:),80,c,'filled'); hold on; scatter(ROI_ax,ROI_values(2,:),100,c,'filled','d');
coeff = corrcoef([ROI_ax ROI_ax], [ROI_values(1,:) ROI_values(2,:)]);
legend(['r = ' num2str(round(coeff(1,2),2))],'location','NorthWest');


set(gca,'FontSize',22)
ylabel('TDR');
xlabel('Axon diameter (\mum)');
title('G_{max} = 2500 mT/m')

%% plot of the ROI definition and TDR map


% pick which of the gradient settings to choose for plotting;
ids = {'Aeon_SC_New_1_G600','Aeon_SC_New_2_G600'};

%ids = {'Aeon_SC_New_1_Gmax','Aeon_SC_New_2_Gmax'};
 
i = 2;
    
mouse_id = ids{i};
data_folder = [filepath '/data/' mouse_id];
load([data_folder filesep mouse_id '_TDR']);
load([data_folder filesep mouse_id '_Den']);
load([data_folder filesep mouse_id '_mask_WM_ROIs']);
load([data_folder filesep mouse_id '_mask_GM']);

if contains(mouse_id,'_2')
    slice = 5;
else 
    slice = 3;
end


figure();
  imagesc(image_den(:,:,slice,1).*mask(:,:,slice)); colormap(gray); axis off; axis square; hold on;
  for j = 1:6
  contour(squeeze(mask_WM_ROI(:,:,slice,j)),'Linecolor',c(j,:),'Linewidth',1); axis square; hold on;
  end

  


figure();
TDR_plot = TDR_opt(:,:,slice,12);
TDR_plot(mask(:,:,slice) == 0) = -2;
cmap = parula(256);
cmap(1,:) = 0;
imagesc(TDR_plot,[-0.1 0.6]); axis square; axis off; colormap(cmap); colorbar;
set(gca,'FontSize',18)
