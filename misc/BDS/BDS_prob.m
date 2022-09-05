function  final_streamlines=BDS_prob(atlas_base,conn_map,map,stride_size,seeding,theta_thresh,connection,N_core)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%This function is used for mapping connectome using  probabilistic BDS
%%%%%algorithm

%%%inputs
%atlas_base: cortical parcellation used for connectome mapping
%conn_lookup: lookup table for intra-block connections
%block_loc: lookup table for block-image
%conn_map: block-image
%stride_size: stride used for craeting block-image
%seeding: seeding methodology (ROI or whole-image)
%theta_thresh: angle threshold

%%%outputs
%conn: connectome reconstructed by the algorithm
%track: complete set of block-chains generated by the algorithm
%final_chains: block-chains only connecting pairs of regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
%connecting regions
p1=connection(1); p2=connection(2);
final_streamlines=[];


%generate lookup tables
[conn_lookup,block_loc]=initilization_variables(atlas_base);

%misc variables
[connec_table,block_table,dir] = find(conn_map');
block_size=4; %block size NxNxN
atlas=create_block_atlas_3d(block_size); %create block-level atlas
nodes=max(max(max(atlas)));
index_table=find(triu(ones(nodes,nodes)));
block_threshold=3; %block_chains should have more than the specified blocks

%look up table for nodes
for i=1:nodes
    [atlas_table(i,1),atlas_table(i,2),atlas_table(i,3)]=ind2sub([block_size,block_size,block_size],find(atlas==i)); 
end

%block-image size
blocks_per_x=size(atlas_base,1)-block_size+1;
blocks_per_y=size(atlas_base,2)-block_size+1;
blocks_per_z=size(atlas_base,3)-block_size+1;

%output file
final_chains=[];

tic
percent=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
final=floor(length(dir).*percent);
frst=0; total=0;



parfor i=1:length(dir)
       check=(i==final);
    %if length(find(check)) >0
    %    disp(['Connections traversed: ' num2str(percent(find(check))*100) ' %']);
    %end
    
    
    %indexes of block and connnectivity within a block
    conn_index=index_table(connec_table(i));
    cx=conn_lookup(conn_index,1); 
    cy=conn_lookup(conn_index,2);
    
    ix=block_loc(block_table(i),1); 
    iy=block_loc(block_table(i),2); 
    iz=block_loc(block_table(i),3);
    
    tx=ix; ty=iy; tz=iz;%saving location for back-tracking
    found=0;%false condition to terminate streamline
    st=block_table(i);%initializing block-chain for each seeding connection
    first=1; % aligining vectors for forward/back-tracking
    org_st=map(st,connec_table(i));
   
    %Block positions 
    start_x=((ix-1)*stride_size)+1; start_y=((iy-1)*stride_size)+1;
    start_z=((iz-1)*stride_size)+1;
    
    %voxels of the block for ROI based seeding
    nz = atlas_base(start_x:start_x+block_size-1, start_y:start_y+block_size-1,start_z:start_z+block_size-1);
   
    if (strcmp(seeding,'roi') && ~isempty(find(nz,1))) || strcmp(seeding,'full')
    
    %forward tracking
    while found==0
        
        %for aliging vector in two directions for forward and back-tracking
        if first==1
            
            origin_vector=cal_orientation(cx,cy,atlas_table);
            raw_ix=ix;
            raw_iy=iy; raw_iz=iz; %for finding the neighborhood block
             
            first=0;
        end
        
        %searching for neighborhood block
        raw_ix=origin_vector(1)+raw_ix;
        raw_iy=origin_vector(2)+raw_iy;
        raw_iz=origin_vector(3)+raw_iz;
        
        %location of neighborhood block
        new_ix=round(raw_ix); new_iy=round(raw_iy);
        new_iz=round(raw_iz);
        
        if new_ix>0 && new_ix<=blocks_per_x && new_iy>0 && new_iy<=blocks_per_y && new_iz>0 && new_iz<=blocks_per_z
           if ~(ix==new_ix && iy==new_iy && iz==new_iz)  
               %extracting location and connections in next block
                new_block=sub2ind([blocks_per_x,blocks_per_y,blocks_per_z],new_ix,new_iy,new_iz);
                    if ~isempty(find(conn_map(new_block,:),1))
                        temp_conn=find(conn_map(new_block,:));
                        ang_comp=inf(length(temp_conn),1);
                        temp_vector=zeros(length(temp_conn),3);
                        
                        %vector conversions of connections in next block &
                %angle calculation w.r.t current vector
                            for multi=1:length(temp_conn)
                                
                                temp_cx=conn_lookup(index_table(temp_conn(multi)),1); 
                                temp_cy=conn_lookup(index_table(temp_conn(multi)),2);
                                
                                
                                temp_vector(multi,:)=cal_orientation(temp_cx,temp_cy,atlas_table);
                                 if (temp_vector(multi,:)*origin_vector')<0
                                   temp_vector(multi,:)=-temp_vector(multi,:);
                                 end
                                 ang_comp(multi)=cal_angle(origin_vector,temp_vector(multi,:));
                            end
                          
                            %randomly selecting a connection that satisfies
                            %the specified angle threshold
                        angle_index1 = find(ang_comp<=theta_thresh);
							if ~isempty(angle_index1)
							
                            angle_index2=randi(length(angle_index1));
							new_cx=conn_lookup(index_table(temp_conn(angle_index1(angle_index2))),1); 
							new_cy=conn_lookup(index_table(temp_conn(angle_index1(angle_index2))),2);
							cx=new_cx; sy=new_cy;
							origin_vector=temp_vector(angle_index1(angle_index2),:);
							
							st=[st; new_block];
                            org_st=[org_st; map(new_block,temp_conn(angle_index1(angle_index2)))];
							
                    else
                            found=1; %termination condition (non of the connections satisfy angle threshold)
                 
                        end
                    else
                    found=1; %termination condition (next block is empty-no connections)
                    end
           end
            ix=new_ix; iy=new_iy; iz=new_iz; %updating location of current block 
        else
            found=1;%termination condition (reached end of image)
        end
    end%forward-tracking
    
    %back-tracking (same as forward tracking but with reverse current orientation)
    ix=tx; iy=ty; iz=tz;
    
    cx=conn_lookup(conn_index,1); 
    cy=conn_lookup(conn_index,2);
    
    first=1; found=0;

    
    while found==0
        if first==1
            origin_vector=cal_orientation(cy,cx,atlas_table);
            raw_ix=ix;
            raw_iy=iy; raw_iz=iz;
            first=0;
        end
        
        raw_ix=origin_vector(1)+raw_ix;
        raw_iy=origin_vector(2)+raw_iy;
        raw_iz=origin_vector(3)+raw_iz;
        new_ix=round(raw_ix); new_iy=round(raw_iy); new_iz=round(raw_iz);
         
        if new_ix>0 && new_ix<=blocks_per_x && new_iy>0 && new_iy<=blocks_per_y && new_iz>0 && new_iz<=blocks_per_z
           if ~(ix==new_ix && iy==new_iy && iz==new_iz)  
                new_block=sub2ind([blocks_per_x,blocks_per_y, blocks_per_z],new_ix,new_iy, new_iz);
                if ~isempty(find(conn_map(new_block,:),1))
                    
                    temp_conn=find(conn_map(new_block,:));
                    ang_comp=inf(length(temp_conn),1);
                    temp_vector=zeros(length(temp_conn),3);
                    for multi=1:length(temp_conn)
                        
                        temp_cx=conn_lookup(index_table(temp_conn(multi)),1); 
                        temp_cy=conn_lookup(index_table(temp_conn(multi)),2);
                        
                        temp_vector(multi,:)=cal_orientation(temp_cx,temp_cy,atlas_table);
                        if (temp_vector(multi,:)*origin_vector')<0
                            temp_vector(multi,:)=-temp_vector(multi,:);
                        end
                        ang_comp(multi)=cal_angle(origin_vector,temp_vector(multi,:));
                    end
                    angle_index1 = find(ang_comp<=theta_thresh);
							if ~isempty(angle_index1)
							
                            angle_index2=randi(length(angle_index1));
							new_cx=conn_lookup(index_table(temp_conn(angle_index1(angle_index2))),1); 
							new_cy=conn_lookup(index_table(temp_conn(angle_index1(angle_index2))),2);
							cx=new_cx; sy=new_cy;
							origin_vector=temp_vector(angle_index1(angle_index2),:);
							
							st=[new_block; st];
                            org_st=[map(new_block,temp_conn(angle_index1(angle_index2))); org_st];
							
                    else
                        found=1;
                 
                    end
                else
                    found=1;
                end
           end
           ix=new_ix; iy=new_iy; iz=new_iz;
        else
            found=1;
        end
    end
    
    
    if length(st)>block_threshold
    %start & end block
        start=st(1);
        stop=st(end);
        [start_indx, start_indy,start_indz]=ind2sub([ blocks_per_x, blocks_per_y, blocks_per_z],start);
        start_x=((start_indx-1)*stride_size)+1; start_y=((start_indy-1)*stride_size)+1;
        start_z=((start_indz-1)*stride_size)+1;
    
        nz = atlas_base(start_x:start_x+block_size-1, start_y:start_y+block_size-1, start_z:start_z+block_size-1);
        ind_start=mode(nz(nz~=0));

        %mapping blocks to voxels for "stop"  (other end-point)  
        [stop_indx, stop_indy, stop_indz]=ind2sub([ blocks_per_x, blocks_per_y, blocks_per_z],stop);
        stop_x=((stop_indx-1)*stride_size)+1; stop_y=((stop_indy-1)*stride_size)+1;
        stop_z=((stop_indz-1)*stride_size)+1;
        nz=atlas_base(stop_x:stop_x+block_size-1, stop_y:stop_y+block_size-1, stop_z:stop_z+block_size-1);
        ind_stop=mode(nz(nz~=0));
     
   
        if ~isempty(ind_start) && ~isempty(ind_stop) && ind_start~=0 && ...
            ind_stop~=0 && ind_start~=ind_stop && (ind_start==p1 || ind_start==p2) && ...
            (ind_stop==p1 || ind_stop==p2)
 
            %total=total+1;
            temporary=mat2cell(full(org_st)',1,size(full(org_st)',2));
            %final_streamlines{total}=full(org_st);
            final_streamlines=[final_streamlines, temporary];
        end
   
    
    end
    end%if for seeding from nodes block only
    
end

time=toc;
disp(['Block-chains Generated : ' num2str(time) ' sec.']);


catch
    disp('Error in Prob BDS: Check input files');
end

