function voxel_radius = radius_voxel_cluster(atlas,Clusters)
%How to call function
%radius_voxel_cluster(atlas,Clus_init)

%Let us say that a cluster of streamlines is connecting regions R1 and R2

for i = 1:size(Clusters)
    %{
 voxel_radius(i,1)=nthroot(size(find(atlas ==Clusters(i,1)),1)*3/(4*3.14),3);%This is the voxel radius of region R1
 voxel_radius(i,2)=nthroot(size(find(atlas ==Clusters(i,2)),1)*3/(4*3.14),3);%This is the voxel radius of region R2
%}
temp=zeros(size(atlas));
temp(find(atlas ==Clusters(i,1)))=1;
s=regionprops3(temp,"PrincipalAxisLength");
voxel_radius(i,1)=mean(s.PrincipalAxisLength)/2;
temp=zeros(size(atlas));
temp(find(atlas ==Clusters(i,2)))=1;
s=regionprops3(temp,"PrincipalAxisLength");
voxel_radius(i,2)=mean(s.PrincipalAxisLength)/2;
end
end