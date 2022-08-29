%%
clear;
clc;

%% Specification of the locations of material points
length = 1.0;
% length: Total length of the bar
ndivx = 1000;
% ndivx: Number of divisions in x direction - except boundary region
dx = length / ndivx;
% dx: Spacing between material points
nbnd = 3;
% nbnd: Number of divisions in the boundary region
totnode = ndivx + nbnd;
% totnode: Total number of material points

coord = zeros(totnode, 1);
% coord: Material point locations

nnum = 0;
% nnum: Material point number

% Material points of the bar
for i = 1:ndivx
    coordx = (dx / 2.0) + (i - 1) * dx;
    nnum = nnum + 1;
    coord(nnum, 1) = coordx;
end

% Material points of the constrained region
for i = 1:nbnd
    coordx = (-1.0 / 2.0 * dx) - (i - 1) * dx;
    nnum = nnum + 1;
    coord(nnum, 1) = coordx;
end

%% Determination of material points inside the horizon of each material point
delta = 3.015 * dx;
% delta: Horizon


pointfam = int32(zeros(totnode, 1));
% pointfam: index array to find the family members in nodefam array
numfam = int32(zeros(totnode, 1));
% numfam: Number of family members of each material point
nodefam = int32(zeros(1000000, 1));
% nodefam: array containing family members of all material points

for i = 1:totnode
    if (i == 1)
        pointfam(i, 1) = 1;
    else
        pointfam(i, 1) = pointfam(i - 1, 1) + numfam(i - 1, 1);
    end
    for j = 1:totnode
        idist = norm(coord(j, :)-coord(i, :));
        if (i ~= j)
            if (idist <= delta)
                numfam(i, 1) = numfam(i, 1) + 1;
                nodefam(pointfam(i, 1)+numfam(i, 1)-1, 1) = j;
            end
        end
    end
end

%% Determination of surface correction factors
radij = dx / 2.0;
% radij: Material point radius
area = dx * dx;
% area: Cross-sectional area
vol = area * dx;
% vol: Volume of a material point

dens = 7850.0;
% dens: Density
emod = 200.0e9;
% emod: Elastic modulus
bc = 2.0 * emod / (area * (delta^2));
% bc: Bond constant


disp = zeros(totnode, 1);
% disp: displacement of a material point
stendens = zeros(totnode, 1);
% stendens: strain energy of a material point
fncst = ones(totnode, 1);

% Loading 1
sedload1 = 0.5 * emod * 1.0e-6;
% sedload1: strain energy density of a material point for the first loading condition based on classical continuum mechanics

for i = 1:totnode
    disp(i, 1) = 0.001 * coord(i, 1);
end

for i = 1:totnode
    stendens(i, 1) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        idist = norm(coord(cnode, :)-coord(i, :));
        nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));

        if (idist <= delta - radij)
            fac = 1.0;
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij);
        else
            fac = 0.0;
        end

        stendens(i, 1) = stendens(i, 1) + 0.5 * 0.5 * bc * ((nlength - idist) / idist)^2 * idist * vol * fac;
    end
    % Calculation of surface correction factor in x direction
    % by finding the ratio of the analytical strain energy density value
    % to the strain energy density value obtained from PD Theory
    fncst(i, 1) = sedload1 / stendens(i, 1);
end

%% Initialization of displacements and velocities
vel = zeros(totnode, 1);
disp = zeros(totnode, 1);

%% Stable mass vector computation
dt = 1.0;
% dt: Time interval

massvec = zeros(totnode, 1);
% massvec: massvector for adaptive dynamic relaxation

for i = 1:totnode
    massvec(i, 1) = 0.25 * dt * dt * (2.0 * area * delta) * bc / dx * 5.0;
end

%% Applied loading
appres = 200.0e6;
% appres: Applied pressure

bforce = zeros(totnode, 1);
% bforce: body load acting on a material point

bforce(ndivx, 1) = appres / (dx);

%% Boundary condition - Constrained region
for i = (ndivx + 1):totnode
    vel(i, 1) = 0.0;
    disp(i, 1) = 0.0;
end

%% Time integration
nt = 10000;
% nt: Total number of time step

pforce = zeros(totnode, 1);
% pforce: total peridynamic force acting on a material point
pforceold = zeros(totnode, 1);
% pforce: total peridynamic force acting on a material point
acc = zeros(totnode, 1);
% acc: acceleration of a material point

velhalf = zeros(totnode, 1);
velhalfold = zeros(totnode, 1);
% vel: velocity of a material point

coord_disp_pd_nt = zeros(ndivx, 3);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
center_node = zeros(nt, 2);
% Peridynamic displacement at center node of all time step.
cn = 0.0;
cn1 = 0.0;
cn2 = 0.0;
for tt = 1:nt
    fprintf("%d/%d\n", tt, nt);

    for i = 1:totnode
        pforce(i, 1) = 0.0;
        for j = 1:numfam(i, 1)
            cnode = nodefam(pointfam(i, 1)+j-1, 1);
            idist = norm(coord(cnode, :)-coord(i, :));
            nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));

            % Volume correction
            if (idist <= delta - radij)
                fac = 1.0;
            elseif (idist <= delta + radij)
                fac = (delta + radij - idist) / (2.0 * radij);
            else
                fac = 0.0;
            end

            % Determination of the surface correction between two material points
            scr = (fncst(i, 1) + fncst(cnode, 1)) / 2.0;

            % Calculation of the peridynamic force in x direction
            % acting on a material point i due to a material point j
            dforce1 = bc * (nlength - idist) / idist * vol * scr * fac * (coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1)) / nlength;

            pforce(i, 1) = pforce(i, 1) + dforce1;
        end
    end

    % Adaptive dynamic relaxation ⬇⬇⬇

    for i = 1:ndivx
        if (velhalfold(i, 1) ~= 0.0)
            cn1 = cn1 - disp(i, 1) * disp(i, 1) * (pforce(i, 1) / massvec(i, 1) - pforceold(i, 1) / massvec(i, 1)) / (dt * velhalfold(i, 1));
        end
        cn2 = cn2 + disp(i, 1) * disp(i, 1);
    end

    if (cn2 ~= 0.0)
        if ((cn1 / cn2) > 0.0)
            cn = 2.0 * sqrt(cn1/cn2);
        else
            cn = 0.0;
        end
    else
        cn = 0.0;
    end

    if (cn > 2.0)
        cn = 1.9;
    end

    for i = 1:ndivx
        % Integrate acceleration over time.
        if (tt == 1)
            velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
        else
            velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
        end

        vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
        disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;

        velhalfold(i, 1) = velhalf(i, 1);
        pforceold(i, 1) = pforce(i, 1);
    end

    % Adaptive dynamic relaxation ⬆⬆⬆

    if (tt == nt)
        for i = 1:ndivx
            coord_disp_pd_nt(i, 1:3) = [coord(i, 1), disp(i, 1), 0.001d0 * coord(i, 1)];
        end
    end

    center_node(tt, 1:2) = [tt, disp(500, 1)];
end

%% plot
subplot(121)
plot(center_node(:, 1), center_node(:, 2));
subplot(122)
hold on;
plot(coord_disp_pd_nt(:, 1), coord_disp_pd_nt(:, 2));
plot(coord_disp_pd_nt(:, 1), coord_disp_pd_nt(:, 3));
legend("Peridynamics", "Analytical");
hold off;