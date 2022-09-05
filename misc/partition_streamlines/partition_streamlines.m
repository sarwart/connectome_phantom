function [new_clus,Clus_init,new_CM, CM, clusters]=partition_streamlines(streamlines,base_atlas,st_thresh,wm, auto, density)



%mask for verifying replication retains within the brain mask
bd_mask=wm+base_atlas; 
bd_mask(find(bd_mask))=1; 
check_atlas=base_atlas;


%find connected streamlines
fprintf('Generating Initial Connectivity Matrix \n')
tic;
[connected_st,connected_rg,Clus_init,CM]=connectivity_matrix(streamlines,base_atlas,auto,st_thresh);
time=toc;
fprintf('Initial Connectivity Matrix Generated:  %f\n',time);



tic
%This is for computing the middle/central streamline in a cluster 
[~,count_init,clusters] = find_streamline_centre(streamlines,Clus_init,connected_rg,connected_st);



%adjusting the connection density
dense_conn=max(max(CM))*0.5;%shortlisting the connections to be adjusted
dense=randi(density,1,length(count_init));
norm_ratio=zeros(size(connected_st));

for c=1:length(count_init)
    if count_init(c)<=dense_conn

        temp_ind=find(connected_rg(:,1)==Clus_init(c,1) & connected_rg(:,2)==Clus_init(c,2));
        norm_ratio(temp_ind)=dense(c);
        temp_ind=find(connected_rg(:,2)==Clus_init(c,1) & connected_rg(:,1)==Clus_init(c,2));
        norm_ratio(temp_ind)=dense(c);
    end
    
    

end

% sanity cheeck whether connection exists after thresholding and adjustment

if  exist('clusters','var')
 
%original connected
for i = 1:size(Clus_init,1) %do this for all connections
    


%replicate based on the randomly sampled neighbourhood
step=norm_ratio(i);


ngh=[];
for k=-step:step
    for l=-step:step
        for m=-step:step
        
        if ~(k==0 && l==0 && m==0)
                ngh=[ngh;[k,l,m]]; 
         end
    
        end
    end
end


new_clus{i}=[];
clus=[];
if step~=0
    
    tic
    count_n=0;
     
       
    check_atlas(find(base_atlas))=2;
    check_atlas(find(base_atlas))=2;
    check_atlas(find(base_atlas==Clus_init(i,1)))=1;
    check_atlas(find(base_atlas==Clus_init(i,2)))=1;
    

    
    for j=1:length(clusters{i})
        
   
        for k=1:length(ngh)
        

            
        new_streamline=(clusters{i}{j} + repmat(ngh(k,:),size(clusters{i}{j},1),1));
        ngh_new=floor(new_streamline);
        ind=sub2ind(size(check_atlas),ngh_new(:,1),ngh_new(:,2),ngh_new(:,3));
        bt=ind(2:end-1);
        bd_check=length(find(bd_mask(ind)));
        if (check_atlas(ind(1))==1 && check_atlas(ind(end))==1 && ...
           length(find(check_atlas(bt)==2))==0 && bd_check==length(ind)) 
           
            count_n=count_n+1;
            clus{count_n}=new_streamline;
            
        end
        end
        
    end

end

new_clus{i}=[clus,clusters{i}]; 


time=toc;
fprintf('Streamlines replicated for Original Connections %d : %f\n',i, time);
end

tic

%create a new connectivity matrix after adjustment
new_CM=zeros(max(max(max(base_atlas))),max(max(max(base_atlas))));

for i=1:size(Clus_init,1)
   
        new_CM(Clus_init(i,1),Clus_init(i,2))=length(new_clus{i});
        new_CM(Clus_init(i,2),Clus_init(i,1))=length(new_clus{i});

end


time=toc;
fprintf('Extracting connections that satisfy the streamline threshold: %f\n', time);

else
    disp('NO CONNECTING STRAEMLINES FOUND: CHECK INPUT SETTINGS');
    new_clus =[]; Clus_init=[];CM=[];new_CM=[];  clusters=[];

end

