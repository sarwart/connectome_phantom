function vector=cal_orientation(c1,c2,atlas_table)
%%%%calculate a unit vector from a inta-block connection
vector=(atlas_table(c1,:)-atlas_table(c2,:))/norm(atlas_table(c1,:)-atlas_table(c2,:)); 
end
