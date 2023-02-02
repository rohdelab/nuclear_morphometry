function coor = fromInd2Coord(ind,Ny)

coor(1,:) = fix(ind/Ny)+1';
coor(2,:) = mod(ind,Ny)';



