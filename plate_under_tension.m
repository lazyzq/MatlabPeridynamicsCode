%%
clear;
clc;

%% Specification of the locations of material points
length = 1.0;
% length: Total length of the bar
width = 0.5;
% width: Total width of the plate
ndivx = 100;
% ndivx: Number of divisions in x direction - except boundary region
ndivy = 50;
% ndivy: Number of divisions in y direction - except boundary region
dx = length / ndivx;
% dx: Spacing between material points
nbnd = 0;
% nbnd: Number of divisions in the boundary region
totnode = ndivx * ndivy + nbnd;
% totnode: Total number of material points

coord = zeros(totnode, 2);
% coord: Material point locations

nnum = 0;
% nnum: Material point number

% Material points of the plate
for i = 1:ndivx
    for j = 1:ndivy
        coordx = -1.0 / 2.0 * length + (dx / 2.0) + (i - 1) * dx;
        coordy = -1.0 / 2.0 * width + (dx / 2.0) + (j - 1) * dx;
        nnum = nnum + 1;
        coord(nnum, 1) = coordx;
        coord(nnum, 2) = coordy;
    end
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
        %         idist = norm(coord(j, :)-coord(i, :));
        idist = sqrt((coord(j, 1) - coord(i, 1))^2+(coord(j, 2) - coord(i, 2))^2);
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
thick = dx;
% thick: Total thickness of the plate

dens = 7850.0;
% dens: Density
emod = 200.0e9;
% emod: Elastic modulus
bc = 9.0 * emod / (pi * thick * (delta^3));
% bc: Bond constant


disp = zeros(totnode, 2);
% disp: displacement of a material point
stendens = zeros(totnode, 2);
% stendens: strain energy of a material point
fncst = ones(totnode, 2);
% fncst: surface correction factors of a material point, 1:loading 1, 2:loading 2

% Loading 1
sedload1 = 9.0 / 16.0 * emod * 1.0e-6;
% sedload1: strain energy density of a material point for the first loading condition based on classical continuum mechanics

for i = 1:totnode
    disp(i, 1) = 0.001 * coord(i, 1);
    disp(i, 2) = 0.0;
end

for i = 1:totnode
    stendens(i, 1) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        %         idist = norm(coord(cnode, :)-coord(i, :));
        %         nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2);

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

% Loading 2
sedload2 = 9.0 / 16.0 * emod * 1.0e-6;
% sedload2: strain energy density of a material point for the second loading condition based on classical continuum mechanics

for i = 1:totnode
    disp(i, 1) = 0.0;
    disp(i, 2) = 0.001 * coord(i, 2);
end

for i = 1:totnode
    stendens(i, 2) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        %         idist = norm(coord(cnode, :)-coord(i, :));
        %         nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2);

        if (idist <= delta - radij)
            fac = 1.0;
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij);
        else
            fac = 0.0;
        end

        stendens(i, 2) = stendens(i, 2) + 0.5 * 0.5 * bc * ((nlength - idist) / idist)^2 * idist * vol * fac;
    end
    % Calculation of surface correction factor in x direction
    % by finding the ratio of the analytical strain energy density value
    % to the strain energy density value obtained from PD Theory
    fncst(i, 2) = sedload2 / stendens(i, 2);
end

%% Initialization of displacements and velocities
vel = zeros(totnode, 2);
disp = zeros(totnode, 2);

%% Stable mass vector computation
dt = 1.0;
% dt: Time interval

massvec = zeros(totnode, 2);
% massvec: massvector for adaptive dynamic relaxation

for i = 1:totnode
    % 5 is a safety factor
    massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc / dx * 5.0;
    massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc / dx * 5.0;
end

%% Applied loading - Left and Right
appres = 200.0e6;
% appres: Applied pressure

bforce = zeros(totnode, 2);
% bforce: body load acting on a material point

% left
for i = 1:ndivy
    bforce(i, 1) = -1.0 * appres / (dx);
end

% right
for i = (totnode - ndivy + 1):totnode
    bforce(i, 1) = appres / (dx);
end

%% Time integration
nt = 1000;
% nt: Total number of time step
alpha = 23.0e-6;
% alpha: Coefficient of thermal expansion
dtemp = 0.0;
% dtemp: Temperature change
pratio = 1.0 / 3.0;
% pratio: Poisson's ratio


pforce = zeros(totnode, 2);
% pforce: total peridynamic force acting on a material point
pforceold = zeros(totnode, 2);
% pforceold: total peridynamic force acting on a material point in the previous time step 1:x-coord, 2:y-coord
acc = zeros(totnode, 2);
% acc: acceleration of a material point

velhalf = zeros(totnode, 2);
velhalfold = zeros(totnode, 2);
% vel: velocity of a material point

coord_disp_pd_nt = zeros(ndivx*ndivy, 4);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
% center_node = zeros(nt, 2);
% % Peridynamic displacement at center node of all time step.
horiCnt = 0;
horizontal_disps = zeros(ndivx, 6);
vertiCnt = 0;
vertical_disps = zeros(ndivy, 6);
steady_check = zeros(nt, 3);

cn = 0.0;
cn1 = 0.0;
cn2 = 0.0;
for tt = 1:nt
    fprintf("%d/%d\n", tt, nt);

    for i = 1:totnode
        pforce(i, 1) = 0.0;
        pforce(i, 2) = 0.0;
        for j = 1:numfam(i, 1)
            cnode = nodefam(pointfam(i, 1)+j-1, 1);
            %             idist = norm(coord(cnode, :)-coord(i, :));
            %             nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
            idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2);
            nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2);

            % Volume correction
            if (idist <= delta - radij)
                fac = 1.0;
            elseif (idist <= delta + radij)
                fac = (delta + radij - idist) / (2.0 * radij);
            else
                fac = 0.0;
            end

            % Determination of the surface correction between two material points
            if (abs(coord(cnode, 2)-coord(i, 2)) <= 1.0e-10)
                theta = 0.0;
            elseif (abs(coord(cnode, 1)-coord(i, 1)) <= 1.0e-10)
                theta = 90.0 * pi / 180.0;
            else
                theta = atan(abs(coord(cnode, 2)-coord(i, 2))/abs(coord(cnode, 1)-coord(i, 1)));
            end
            scx = (fncst(i, 1) + fncst(cnode, 1)) / 2.0;
            scy = (fncst(i, 2) + fncst(cnode, 2)) / 2.0;
            scr = 1.0 / (((cos(theta))^2.0 / (scx)^2.0) + ((sin(theta))^2.0 / (scy)^2.0));
            scr = sqrt(scr);

            % Calculation of the peridynamic force in x direction
            % acting on a material point i due to a material point j
            dforce1 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1)) / nlength;
            dforce2 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2)) / nlength;

            pforce(i, 1) = pforce(i, 1) + dforce1;
            pforce(i, 2) = pforce(i, 2) + dforce2;
        end
    end

    % Adaptive dynamic relaxation ⬇⬇⬇

    for i = 1:totnode
        if (velhalfold(i, 1) ~= 0.0)
            cn1 = cn1 - disp(i, 1) * disp(i, 1) * (pforce(i, 1) / massvec(i, 1) - pforceold(i, 1) / massvec(i, 1)) / (dt * velhalfold(i, 1));
        end
        if (velhalfold(i, 2) ~= 0.0)
            cn1 = cn1 - disp(i, 2) * disp(i, 2) * (pforce(i, 2) / massvec(i, 2) - pforceold(i, 2) / massvec(i, 2)) / (dt * velhalfold(i, 2));
        end
        cn2 = cn2 + disp(i, 1) * disp(i, 1);
        cn2 = cn2 + disp(i, 2) * disp(i, 2);
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

    for i = 1:totnode
        % Integrate acceleration over time.
        if (tt == 1)
            velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
            velhalf(i, 2) = 1.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2)) / 2.0;
        else
            velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
            velhalf(i, 2) = ((2.0 - cn * dt) * velhalfold(i, 2) + 2.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2))) / (2.0 + cn * dt);
        end

        vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
        vel(i, 2) = 0.5 * (velhalfold(i, 2) + velhalf(i, 2));
        disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;
        disp(i, 2) = disp(i, 2) + velhalf(i, 2) * dt;

        velhalfold(i, 1) = velhalf(i, 1);
        velhalfold(i, 2) = velhalf(i, 2);
        pforceold(i, 1) = pforce(i, 1);
        pforceold(i, 2) = pforce(i, 2);
    end

    % Adaptive dynamic relaxation ⬆⬆⬆

    if (tt == nt)
        for i = 1:(ndivx * ndivy)
            coord_disp_pd_nt(i, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
            if (abs(coord(i, 2)-(dx / 2.0)) <= 1.0e-8)
                horiCnt = horiCnt + 1;
                horizontal_disps(horiCnt, 1:6) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), 0.001 * coord(i, 1), -1.0 * 0.001 * pratio * coord(i, 2)];
            end
            if (abs(coord(i, 1)-(dx / 2.0)) <= 1.0e-8)
                vertiCnt = vertiCnt + 1;
                vertical_disps(vertiCnt, 1:6) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), 0.001 * coord(i, 1), -1.0 * 0.001 * pratio * coord(i, 2)];
            end
        end
    end

    steady_check(tt, 1:3) = [tt, disp(3788, 1), disp(3788, 2)];
end

%% plot
colormap jet;
subplot(221)
scale = 500;
scatter(coord_disp_pd_nt(:, 1)+scale*coord_disp_pd_nt(:, 3), coord_disp_pd_nt(:, 2)+scale*coord_disp_pd_nt(:, 4), [], sqrt(coord_disp_pd_nt(:, 3).^2+coord_disp_pd_nt(:, 4).^2), "filled")
subplot(222)
plot(steady_check(:, 1), steady_check(:, 2), steady_check(:, 1), steady_check(:, 3))
subplot(223)
plot(horizontal_disps(:, 1), horizontal_disps(:, 3), horizontal_disps(:, 1), horizontal_disps(:, 5))
subplot(224)
plot(vertical_disps(:, 2), vertical_disps(:, 4), vertical_disps(:, 2), vertical_disps(:, 6))
