function g = domainSimpleY()

coordV = ...
    [0, 0; ...
    0, 4; ...
    4, 4; ...
    4, 0; ...
    2, 1; ...
    1, 2; ...
    2, 3; ...
    3, 2];

V0T = ...
    [1, 4, 5; ...
    1, 5, 6; ...
    1, 6, 2; ...
    6, 7, 2; ...
    2, 7, 3; ...
    7, 8, 3; ...
    8, 4, 3; ...
    5, 4, 8];

g = Grid(coordV, V0T);

g.idE(1:4) = 1;
g.idE([9, 12, 15, 16]) = 5;
printline(1, 'Edge IDs have been defined:')
printline(2, 'edges on the interiour boundary: 1')
printline(2, 'edges on the interiour boundary: 5')

end