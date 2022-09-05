function theta=cal_angle(v1,v2)
%%%calculate an angle between two vectors
theta=real(acosd(((v1(1)*v2(1))+(v1(2)*v2(2))+(v1(3)*v2(3)))/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2))));
end              

                