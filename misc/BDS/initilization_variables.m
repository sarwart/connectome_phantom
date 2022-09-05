function [conn_lookup,block_loc]=initilization_variables(atlas_base)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%This function is used to create and initialize lookup tables used by
%%%%%variants of BDS algorithms for mapping connectome 

%%%input
%atlas_base: cortical parcellation used for extracting dimensions

%%%ouputs
%conn_lookup: lookup table for intra-block connections
%block_loc: lookup table for block-image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Creating and Initializing Variables....');

%block-level atlas
block_size=4; %block size NxNxN
atlas=create_block_atlas_3d(block_size); %create block-level atlas
nodes=max(max(max(atlas)));


%block-image size
blocks_per_x=size(atlas_base,1)-block_size+1;
blocks_per_y=size(atlas_base,2)-block_size+1;
blocks_per_z=size(atlas_base,3)-block_size+1;

%%%%Creating Look-up Tables(saves computational time in every iteration)%%%%

tic

%indexes of block-level connectivity matrix into subscripts 
conn_lookup=zeros(nodes^2,2);
for i=1:nodes^2
[conn_lookup(i,1),conn_lookup(i,2)]=ind2sub([nodes,nodes],i);
end

%look-up table for blocks location 

block_loc=zeros(blocks_per_x*blocks_per_y*blocks_per_z,3);

for i=1:(blocks_per_z*blocks_per_y*blocks_per_x)
   
    [block_loc(i,1),block_loc(i,2),block_loc(i,3)]=ind2sub([blocks_per_x,blocks_per_y,blocks_per_z],i);   
    
end

time=toc;
disp(['Look-up Tables Creation: ' num2str(time) ' sec.']);