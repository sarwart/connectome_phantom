function atlas=create_block_atlas_3d(N)
%%%% create a block-level atlas
atlas=ones(N,N,N);
atlas(2:end-1,2:end-1,2:end-1) = 0;
non_zero=find(atlas);
val=1:length(non_zero);
atlas(non_zero)=val;
end