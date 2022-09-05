function [dist,angle] = init_dis_ang_cal(discard,D_n,A_n)


%Using the Avg/midpoint concept
centre = sum(discard,1)/length(discard); %Adding all the points in the discarded streamline and then dividing by length to find center.
%this (N+1)/2 point
dist(1,:) = discard(1,:); %This is the the first point in the discarded streamline
dist(D_n,:) = discard(end,:); %This is the the last point in the discarded streamline

dist((D_n + 1)/2,:) = centre; %This is the the centre point in the discarded streamline

%Computing other points in between first and centre point of the streamline
for i = 2 : (D_n - 1)/2
    dist(i,:) = (i-1)*(dist(1,:) + centre)/((D_n - 1)/2);
end
%Computing other points in between last and centre point of the streamline
for i = (D_n + 3)/2 : (D_n - 1)
    dist(i,:) = (i-(D_n + 1)/2)*(dist(D_n,:) + centre)/((D_n - 1)/2);
end
%Now the matrix dist contains D_n number of points which are equidistant to
%each other

angle(1) = cal_angle_3d(discard(1,:),discard(round(length(discard)/(A_n)),:)); %This is the the first angle in the discarded streamline
angle(A_n) = cal_angle_3d(discard(round((A_n-1)*length(discard)/(A_n)),:),discard(end,:)); %This is the the last angle in the discarded streamline

%Computing other angles present in the discarded streamline the streamline
for i = 2: (A_n-1)
    %Here we use the cal_angle_3d function to find angle between two
    %adjacent points in the discarded streamline
    angle(i) = cal_angle_3d(discard(round((i-1)*length(discard)/(A_n)),:),discard(round((i+1)*length(discard)/(A_n)),:));
end
%Now the matrix angle contains A_n number of angles between consecutive
%points in the discarded streamline
end