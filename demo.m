
clc;
clear;
close all;
addpath(genpath(strcat(pwd,'/misc')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initializing file paths %%%%%%%%%%%%%%%%%%%%%%%

save_path=strcat(pwd,'/data/');

    
%tic
%read path to read tracks

file=strcat(pwd,'/data/streamlines.tck');
%file='E:\Rishabh\21\less_new1_ntracks_21.tck';

%read the cortical parcellation
temp=load_untouch_nii(strcat(pwd,'/data/atlas.nii.gz'));
atlas=temp.img;

%read whole-brain mask
temp=load_untouch_nii(strcat(pwd,'/data/brain_mask.nii.gz'));
nbm=double(temp.img);

%read wm-mask for replication check
temp=load_untouch_nii(strcat(pwd,'/data/WM.nii.gz'));
boundry_mask=double(temp.img);

%read CSF-mask
temp=load_untouch_nii(strcat(pwd,'/data/CSF.nii.gz'));  
csf=double(temp.img);


%read GM-mask
temp=load_untouch_nii(strcat(pwd,'/data/GM.nii.gz'));
gm=double(temp.img);


%gradient file
parameters.g=load(strcat(pwd,'/data/bvecs'));

%single shell with one b-value
constant_b=2000;
parameters.b=ones(size(parameters.g,2),1)*constant_b;
parameters.b(1)=0;


%%%%%%%%%%%%%%%%%%%%%%%% Initializing parameters  %%%%%%%%%%%%%%%%%%%%%%%

%%%dMRI acquisition parameters
output_voxel=1.25; %voxel size for output dMRI data


snr=20; %quality for output dMRI 

%diffusivity for WM (zeppelin model)
L1=2.2*10^-3; 
L2=0.2*10^-3; 

parameters.diff_gm=0.9*10^-3; %diffusivity for GM
parameters.diff_csf=3*10^-3; %diffusivity for CSF

parameters.wm_T2=44*1000; %T2 for WM 
parameters.gm_T2=51*1000; %T2 for GM 
parameters.csf_T2=500*1000; %T2 for CSF 

parameters.wm_T1=832*1000; %T1 for WM 
parameters.gm_T1=1331*1000; %T1 for GM 
parameters.csf_T1=3700*1000; %T1 for CSF 

parameters.TE=57*1000; %TE
parameters.TR=8800*1000; %TR


%these are used only when a voxel also contains WM
parameters.fu_gm=0.35; %volume fraction for GM
parameters.fu_csf=0.35; % volume fraction for CSF


density=[2, 6]; %controlling the fiber density in voxels


auto=0; % if 0 then use manual threshold, if 1 then ut will threshold 5% weak connections
st_thresh=1; %threshold for connecting streamlines (if automated threshold is not used)

%for signal simulation
[r,c,z]=size(atlas);
dwi_empty=zeros(r,c,z,size(parameters.g,2));
dwi=dwi_empty;
mask=zeros(r,c,z);%WM mask for new dMRI

%%%%%%%%%%%%%%%%%%%%%%% Reading tracks/Streamlines %%%%%%%%%%%%%%%%%%%%%%%


%for reading tck file

tracks=read_mrtrix_tracks(file);
streamlines=cell(1,length(tracks.data));

for i=1:length(tracks.data)
streamlines(i)=tracks.data(i);
end
clear tracks.data;
%{

%for reading trk file
[header tracks]=trk_read(file);
streamlines=cell(1,length(tracks));
for i=1:length(tracks)
streamlines{i}=tracks(i).matrix;

end
clear tracks;
header.n_properties=0;
%}


%%%%%%%%%%%%%%%%%%%%%%% Extract streamlines for each connection %%%%%%%%%%%%%%%%%%%%%%%

[new_clus,Clus_init,new_CM, CM, org_clus]=partition_streamlines(streamlines,atlas,st_thresh,boundry_mask, auto, density,1); %step-size of 1 used for reducing the computational time



%%%%%%%%%%%%%%%%%%%%%%% Simulate Signal %%%%%%%%%%%%%%%%%%%%%%%

dense=zeros(r,c,z,size(Clus_init,1));

disp('Simulating signal for each bundle');
for C=1:size(Clus_init,1)
    tic;
    [mask_temp,dwi_temp,dense(:,:,:,C)]=track2dwi(new_clus{C},[r,c,z],parameters.g,parameters.b,L1,L2);
    dwi=dwi+(dwi_temp.*dense(:,:,:,C));
    mask=mask+mask_temp;
    time=toc;
    clear mask_temp dwi_temp;
    fprintf('Signal simulated for Bundle %d : %f\n',C, time);
end

dense=sum(dense,4);
dense(find(dense==0))=1;
dwi=dwi./dense;
tdense=max(max(max(dense)));
dwi=dwi.*(dense/tdense);

    
%to simulate 
disp('Simulating signal for the whole brain');
tic
%computing masks for new dMRI
div=mask;
div(find(div==0))=1;
mask(find(mask))=1;



gm_temp=nbm-mask-csf;
gm_temp(find(gm_temp<0))=0;
gm_temp=gm+gm_temp;
gm_temp(find(gm_temp))=1;

grad0=find(parameters.b~=0);


dwi=simulate_signal(dwi,mask,gm_temp,csf,parameters);


%rextracting mask for computing noise
mask=mask.*boundry_mask;
%gm=gm_temp - mask;

dwi=(1000*output_voxel^3).*dwi;
S0=squeeze(dwi(:,:,:,1)); 

%adding noise
n=mean(S0(find(mask)))/snr; 
dwi=ricernd(dwi,n);

nbm=~nbm;


nbm=repmat(nbm, 1,1,1,size(parameters.g,2));

index=find(nbm);
dwi(index)=0;
nii=make_nii(dwi,output_voxel);
save_nii(nii, strcat(save_path,'dMRI.nii.gz'));

nii=make_nii(mask,output_voxel);
save_nii(nii, strcat(save_path,'new_WM.nii.gz'));

save(strcat(save_path,'variables'),'CM','new_CM','new_clus');

nii=make_nii(atlas,output_voxel);
save_nii(nii,strcat(save_path,'sim_atlas.nii.gz'));


dlmwrite(strcat(save_path,'bvecs'),parameters.g,'delimiter',' ');
dlmwrite(strcat(save_path,'bvals'),parameters.b,'delimiter',' ');

time=toc;

fprintf('Signal simulated for whole-brain %f\n', time);



