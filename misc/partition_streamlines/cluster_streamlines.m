function [st_check,newclus]=cluster_streamlines(streamlines,central,count_init,cluster,Clus_init,discard_check,st_check,angle_clus,voxel_radius,thresh_ang,D_n,A_n)
tic;
for i = 1:size(Clus_init,1)
newclus{i} = cluster{i}; %Initializing the newclus with the initial clusters that we got from connectivity matrix. Later in this matrix, we will add the new clusters that we get after the clustering algorithm
countnew(i) = count_init(i); %countnew matrix keeps the track of number of streamlines in a particular cluster. We first initialize it with the initail count values that we get using the connectivity matrix
end

disp('clustering');
    for g = 1: length(discard_check)
        %To find distance and angles of each discarded streamline
        [dist_disc, angle_disc] = init_dis_ang_cal(streamlines{discard_check(g)},D_n,A_n);
        for C = 1:size(Clus_init)
            
            for p = 1:A_n
                %Check if it satisfies angle condition
                if (abs(angle_clus(C,p)- angle_disc(p)) <=thresh_ang)
                else break %Break the loop if any one angle doesn't satisfy the condition 
                end
            end
            if (p==A_n) %This condition is satisfied if all the A_n satisfy the angle condition. If yes, then move towards next step (checking distance condition)

             
             %Below function finds number of points satisfying the distance
             threshold = threshold_dist(central{C},dist_disc,D_n,voxel_radius(C,:));
             
                if  (threshold >= (D_n + 1)/2)
                   countnew(C) = countnew(C) + 1;
                
                   newclus{C}{countnew(C)} = streamlines{discard_check(g)};
                   st_check(discard_check(g)) = 1; %This matrix keeps track of the streamlines which are connecting a region
                end
            end
        end
    end
