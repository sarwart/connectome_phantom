function [wm,dwi,dense]=track2dwi(tracks,dim,g,b,L1,L2)



wm=zeros(dim);
dwi=zeros([dim, length(g)]);
%read the streamlines and generate signal

for i=1:length(tracks)
    
   
    streamline=tracks{i};
    coord=abs(ceil(tracks{i}));
    
    
    for j=1:length(coord)-1
        u=streamline(j+1,:)-streamline(j,:);
        v=u/norm(u);
        
        D=v'*v*(L1-L2)+eye(3)*L2; 
        %stick and zepplin model
        s=(0.8.*exp(diag(-b*L1*(v*g).^2)))+(0.2.*exp(diag(-b.*(g'*D*g)))); 
    
        if [coord(j,1),coord(j,2),coord(j,3)]<=dim & [coord(j,1),coord(j,2),coord(j,3)]>=[1,1,1]

        dwi(coord(j,1),coord(j,2),coord(j,3),:)=squeeze(dwi(coord(j,1),coord(j,2),coord(j,3),:))+ s;
        wm(coord(j,1),coord(j,2),coord(j,3))=wm(coord(j,1),coord(j,2),coord(j,3))+1;
        
        end

    end

end

dense=wm;
div=wm;
div(find(div==0))=1;
dwi=dwi./div;

wm(find(wm))=1;

