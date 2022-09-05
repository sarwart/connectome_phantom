function tracks=blocks_to_streamlines(atlas_base,final_chains,vox,file_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%This function is used mapping block-chains to streamlines for
%%%%%visualization purposes

%%%inputs
%atlas_base: cortical parcellation used for extracting dimensions
%final_chains: block-chains which needs to be converted into streamlines
%vox: voxel size for mapping  streamlines (corresponding to respective dMRI dataset)
%file_type: output file type (.trk or .tck)

%%%outputs
%tracks: streamlines representation of the input block-chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
%block-image size
block_size=4;
blocks_per_x=size(atlas_base,1)-block_size+1;
blocks_per_y=size(atlas_base,2)-block_size+1;
blocks_per_z=size(atlas_base,3)-block_size+1;


if strcmp(file_type,'.tck')
    %%%creating streamlines that stasify .tck format
    tracks= cell(1,length(final_chains));
for i=1:length(final_chains)
    temp=zeros(length(final_chains{i}),3);
    for j=1:length(final_chains{i})
        [temp(j,1),temp(j,2),temp(j,3)]=ind2sub([ blocks_per_x, blocks_per_y, blocks_per_z],final_chains{i}(j));
    end
    tracks{i}=temp.*vox;
end
elseif strcmp(file_type,'.trk')
    %%%creating streamlines that stasify .trk format
    field1 = 'matrix';  value1 = cell(1, length(final_chains));
    field2 = 'nPoints';  value2 = 0;
    tracks = struct(field1,value1,field2,value2);
    for i=1:length(final_chains)
        temp=zeros(length(final_chains{i}),3);
    for j=1:length(final_chains{i})
        [temp(j,1),temp(j,2),temp(j,3)]=ind2sub([ blocks_per_x, blocks_per_y, blocks_per_z],final_chains{i}(j));
    end
    tracks(i).matrix=temp.*vox;
    tracks(i).nPoints=length(temp);
    end
else
    disp('unrecognized file type');
    tracks=[];
end
catch
    disp('Error in Mapping Blocks to Streamlines: Check input files');
    tracks=[];
end