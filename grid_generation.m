grids.x_u = x_min:dx:x_max;
grids.y_u = y_min+dy/2:dy:y_max-dy/2;

grids.x_v = x_min+dx/2:dx:x_max-dx/2;
grids.y_v = y_min:dy:y_max;

grids.x_p = x_min+dx/2:dx:x_max-dx/2;
grids.y_p = y_min+dy/2:dy:y_max-dy/2;

N.y_u = length(grids.y_u); 
N.x_u = length(grids.x_u);
N.y_v = length(grids.y_v); 
N.x_v = length(grids.x_v);
N.y_p = length(grids.y_p); 
N.x_p = length(grids.x_p);