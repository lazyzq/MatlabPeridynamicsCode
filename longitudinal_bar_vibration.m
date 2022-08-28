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

%% Initial condition
for i = 1:ndivx
    vel(i, 1) = 0.0;
    disp(i, 1) = 0.001 * coord(i, 1);
end

%% Boundary condition - Zero displacement at x = 0
for i = (ndivx + 1):totnode
    vel(i, 1) = 0.0;
    disp(i, 1) = 0.0;
end

%% Time integration
nt = 26000;
% nt: Total number of time step
dt = 0.8 * sqrt(2.0*dens*dx/(2.0 * delta * area * bc));
% dt: Time step size
ntotrao = 20;
% ntotrao: Number of terms in the summation for the analytical displacement calculation
cwave = sqrt(emod/dens);
% cwave: Wave speed

pforce = zeros(totnode, 1);
% pforce: total peridynamic force acting on a material point
bforce = zeros(totnode, 1);
% bforce: body load acting on a material point
acc = zeros(totnode, 1);


andisp = zeros(nt, 1);
% andisp: analytical displacements for results
pddisp = zeros(nt, 1);
% pddisp: peridynamic displacements for results
pdtime = zeros(nt, 1);
% pdtime: time array for results

for tt = 1:nt
    fprintf("%d/%d\n", tt, nt);
    ctime = tt * dt;

    for i = 1:ndivx
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

    for i = 1:ndivx
        acc(i, 1) = (pforce(i, 1) + bforce(i, 1)) / dens;
        % Calculate the acceleration of material point i
        vel(i, 1) = vel(i, 1) + acc(i, 1) * dt;
        % Calculate the displacement of material point i
        % by integrating the velocity of material point i
        disp(i, 1) = disp(i, 1) + vel(i, 1) * dt;
    end

    % Store the displacement and time information for the material point at the center
    % of the bar for results
    pddisp(tt, 1) = disp(500, 1);
    pdtime(tt, 1) = ctime;
    % Calculate the analytical displacement solution of the material point at the center
    % of the bar
    for nrao = 0:ntotrao
        andisp(tt, 1) = andisp(tt, 1) + ((-1.0)^(nrao)) / ((2.0 * nrao + 1.0)^2) * sin((2.0 * nrao + 1.0)*pi*coord(500, 1)/2.0) * cos((2.0 * nrao + 1.0)*pi*cwave*ctime/2.0);
    end
    andisp(tt, 1) = 8.0 * 0.001 * 1.0 / (pi^2) * andisp(tt, 1);
end

%% plot
plot(pdtime(1:nt, 1), pddisp(1:nt, 1), pdtime(1:nt, 1), andisp(1:nt, 1));
legend("Peridynamics", "Analytical")
xlabel("Tims(s)");
ylabel("Displacement(m)")