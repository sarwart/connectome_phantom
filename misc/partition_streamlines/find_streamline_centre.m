function [central,count,cluster] = find_streamline_centre(streamlines,Clusters,conn_region,conn_st)
    
count = zeros(size(Clusters,1),1); %Initializing the matrix with zeros
%Count is a matrix which contains number of streamlines for all clusters

for i = 1:size(Clusters,1) %do this for all clusters
    for b = 1:size(conn_st,2) %do this for all connected streamlines
        if (conn_region(b,:) == Clusters(i,:)) % Here we have this condition because we want to partition clusters which are connecting two particular regions
         count(i) = count(i)+1; 
         cluster{i}{count(i)} = streamlines{conn_st(b)}; %This matrix contains streamlines present in a particular cluster
         NPOINT(count(i)) = length(cluster{i}{count(i)}); %This matrix contains length of each streamline in a cluster
        end
    end
    maxNPOINT(i) = max(NPOINT); %This is equal to the maximum of the length of streamlines present in a cluster
    %We need this matrix because we want to interpolate all the streamlines
    clear NPOINT
    
    central{i} = zeros(maxNPOINT(i),3); %Initialize the central streamline
    
    for b = 1:count(i) %do this for all streamlines present in a cluster 
    
    for j = 1:3 %We interpolate along all axes x,y and z. Thats why j =1,2,3
    clusterinterpol{i}{b}(:,j) = interp1(1:length(cluster{i}{b}(:,j)),cluster{i}{b}(:,j),linspace(1, length(cluster{i}{b}(:,j)), maxNPOINT(i)), 'linear'); %this function interpolates the streamlines
    end
    central{i} = central{i} + clusterinterpol{i}{b}; %We add all the interpolated streamlines
    end
    central{i} = central{i}/count(i); %we divide it by total no. of streamlines to get the central streamline

end
end