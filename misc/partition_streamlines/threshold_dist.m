function threshold = threshold_dist(central,dist_dis,D_n,radius)
threshold = 0;
%How to call this function
%dist_threshold(central(C),dist_dis,D_n,voxel_radius(C,:))

for i=1:D_n
dist = (central - dist_dis(i,:)); %Subtracting (a point 'D' from the discarded streamline) from the central streamline
dist = sqrt(sum(dist.^2,2)); %since we have 3 coordinates. We find Eucledian distance
min_dis = min(dist); %Now we find the min value of the matrix dist
%Thus min_dis gives us the minimum distance between point D and the central
%streamline
if(min_dis <= radius(1) && min_dis <= radius(2)) %If this minimum distance is less than or equal to both the radius of the two regions connected by this central straemline
    threshold = threshold + 1; %Then we increment the threshold
end
end

end