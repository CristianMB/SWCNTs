% ........................... Input Data ...................................
n = input('n Chiral Index: ');
m = input('m Chiral Index: ');
l = input('Tube lenght (nm): ');

%Interatomic Spacing C-C
bl = 0.144;

%Lattice Constant
a = 2 * bl * cos(pi / 6);

a1 = [1; 0] * a;
a2 = [cos(pi / 3); sin(pi / 3)] * a;

vector = n * a1 + m * a2;

length_vector = sqrt(vector(1, 1)^2 + vector(2, 1)^2);
cos_t = vector(1, 1) / length_vector;
sin_t = vector(2, 1) / length_vector;
T = [cos_t, sin_t; -sin_t, cos_t];
counter = 1;
check = 1;

y_temp = -bl / 2;


% ...................... Generating graphene plane .........................

graphene_plane = [];
for j = 0:5 * round(l / bl)
    if (check == 3)
        xo = a / 2;
        check = check + 1;
        y_temp = y_temp + bl / 2;
    elseif (check == 4)
        x = a / 2;
        check = 1;
        y_temp = y_temp + bl;
    elseif (check == 2)
        xo = 0;
        check = check + 1;
        y_temp = y_temp + bl;
    else
        xo = 0;
        check = check + 1;
        y_temp = y_temp + bl / 2;
    end
    for i = -5 * round(l / a):5 * round(l / a)
        x = xo + i * a;
        y = y_temp;
        R = [x; y];
        R_temp = T * R;
        if (and(R_temp(1, 1) >= 0, R_temp(1, 1) < 1.01 * length_vector))
            if (and(R_temp(2, 1) >= 0, R_temp(2, 1) < l))
                graphene_plane(counter, 1) = counter;
                graphene_plane(counter, 2) = R_temp(1, 1);
                graphene_plane(counter, 3) = R_temp(2, 1);
                counter = counter + 1;
            end
        end
    end
end

% .............. Rolling the Graphene plane ..........................
Number_atoms = size(graphene_plane, 1);
r = length_vector / 2 / pi;

Atom_data = zeros(Number_atoms, 4);

for i = 1:Number_atoms
    teta = graphene_plane(i, 2) / r;
    Atom_data(i, 1) = i;
    Atom_data(i, 2) = r * sin(teta);
    Atom_data(i, 3) = graphene_plane(i, 3);
    Atom_data(i, 4) = r * cos(teta) + r;
end

% ............. Generating Bond matrix ...............................
Bond_matrix = zeros(Number_atoms);

for i = 1:Number_atoms
    for j = i + 1:Number_atoms
        distance = norm(Atom_data(i, 2:4) - Atom_data(j, 2:4));

        if (distance < 1.1 * bl)
            Bond_matrix(i, j) = 1;
            Bond_matrix(j, i) = 1;
        end
    end
end

% .................. Drawing atoms and bonds .........................
figure;
hold on;
% Set color for atoms
atom_color = [0 0 0];

for i = 1:Number_atoms
    x = Atom_data(i, 2);
    y = Atom_data(i, 3);
    z = Atom_data(i, 4);
    [xx, yy, zz] = sphere();
    rr = 0.08 * a;

    surf(rr * xx + x, rr * yy + y, rr * zz + z, 'EdgeColor', 'none', 'LineStyle', 'none', 'FaceLighting', 'phong', 'FaceColor', atom_color)
    axis equal;
end

% Set color for bonds
bond_color = [0.50 0.000 0.50];


for i = 1:Number_atoms
    for j = i + 1:Number_atoms
        if (Bond_matrix(i, j) == 1)
            x = Atom_data(i, 2);
            y = Atom_data(i, 3);
            z = Atom_data(i, 4);
            xx = Atom_data(j, 2);
            yy = Atom_data(j, 3);
            zz = Atom_data(j, 4);

            line([x, xx], [y, yy], [z, zz], 'Color', bond_color, 'LineWidth', 1)
        end
    end
end

% .................. Generating output files .........................
Atom_data = Atom_data';
Bond_data = find(Bond_matrix == 1);
[i, j] = ind2sub(size(Bond_matrix), Bond_data);
Bond_data = [i, j];

fid1 = fopen('Atom_data.txt', 'wt');
fid2 = fopen('Bond_data.txt', 'wt');
fprintf(fid1, 'Atom_number     X(nm)     Y(nm)     Z(nm)     \n');
fprintf(fid2, 'Bond_number      Atom1      Atom2     \n');
fprintf(fid1, '%d  %15.10f  %15.10f  %15.10f  \n', Atom_data);
fprintf(fid2, '%d  %d  %d  \n', Bond_data);
fclose(fid1);
fclose(fid2);

hold off;
