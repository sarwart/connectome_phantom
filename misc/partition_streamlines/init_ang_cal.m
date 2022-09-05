function angle = init_ang_cal(central,C_n,A_n)
%This function is for finding the angles of the central streamlines of each
%cluster

for C = 1:C_n
angle(C,1) = cal_angle_3d(central{C}(1,:),central{C}(round(length(central{C})/(A_n)),:)); %This is the the first angle
angle(C,A_n) = cal_angle_3d(central{C}(round((A_n-1)*length(central{C})/(A_n)),:),central{C}(end,:)); %This is the the last angle

for i = 2: (A_n-1)
    %Here we use the cal_angle_3d function to find angle between two
    %adjacent/consecutive points in the central streamline
    angle(C,i) = cal_angle_3d(central{C}(round((i-1)*length(central{C})/(A_n)),:),central{C}(round((i+1)*length(central{C})/(A_n)),:));
end
end
end