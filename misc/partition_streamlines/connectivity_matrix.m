function [connected_st,connected_rg,Clus_init,CM]=connectivity_matrix(streamlines,atlas, auto, thresh)

tic
%%compare the tracks with atlas
connected=0;
matrix=zeros( max(max(max(atlas))), max(max(max(atlas))));
CM=matrix;

for t=1:length(streamlines)
    
   % normalizing the streamlines for finding the connectivity matrix
   first=floor(streamlines{t}(1,:));
   last=floor(streamlines{t}(end,:));
   p1=atlas(first(1,1),first(1,2),first(1,3));
   p2=atlas(last(1,1),last(1,2),last(1,3));
   
 %%calculate connectivity matrix
   if p1~=0 && p2~=0 && (p1~=p2)
   matrix(p1,p2)=matrix(p1,p2)+1;
   
   CM(p1,p2)=CM(p1,p2)+1;
   CM(p2,p1)=CM(p2,p1)+1;
   
   connected=connected+1;
   connected_st(connected)=t; %This matrix contains connected streamlines that we get from connectivity matrix
   connected_rg(connected,1) = p1; %This matrix contains the regions that the connected strealines are connecting
   connected_rg(connected,2) = p2;
    end
end

%calculate new threshold if auto threshold is used, else use provided
%thresh
if auto==1
    thresh=floor(0.05 * (max(max(CM))));
end

%threshold on number of connections
CM(find(CM<thresh))=0;

%Finding out number of connections
[row,col] = find(matrix>=thresh);

Clus_init = [row,col];
