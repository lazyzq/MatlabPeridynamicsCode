%%
clear;
clc;

%% Specification of the locations of material points
length = 1.0;
% length: Total length of the block
width = 0.1;
% width: Total width of the block
thick = 0.1;
% thick: Total thickness of the block
ndivx = 100;
% ndivx: Number of divisions in x direction - except boundary region
ndivy = 10;
% ndivy: Number of divisions in y direction - except boundary region
ndivz = 10;
% ndivz: Number of divisions in z direction - except boundary region
dx = length / ndivx;
% dx: Spacing between material points
nbnd = 3;
% nbnd: Number of divisions in the boundary region
totnode = (ndivx + nbnd) * ndivy * ndivz;
% totnode: Total number of material points

coord = zeros(totnode, 3);
% coord: Material point locations
alflag = zeros(totnode, 1);
% flag the points in boundary zone


nnum = 0;
% nnum: Material point number

% Material points of the block
for i = 1:ndivz
    for j = 1:ndivy
        for k = 1:ndivx
            coordx = (dx / 2.0) + (k - 1) * dx;
            coordy = -1.0 / 2.0 * width + (dx / 2.0) + (j - 1) * dx;
            coordz = -1.0 / 2.0 * thick + (dx / 2.0) + (i - 1) * dx;
            nnum = nnum + 1;
            coord(nnum, 1) = coordx;
            coord(nnum, 2) = coordy;
            coord(nnum, 3) = coordz;
            if (coordx > (length - dx))
                alflag(nnum, 1) = 1;
            end
        end
    end
end

totint = nnum;

% Material points of the boundary region - left
for i = 1:ndivz
    for j = 1:ndivy
        for k = 1:nbnd
            coordx = -(dx / 2.0) - (k - 1) * dx;
            coordy = -1.0 / 2.0 * width + (dx / 2.0) + (j - 1) * dx;
            coordz = -1.0 / 2.0 * thick + (dx / 2.0) + (i - 1) * dx;
            nnum = nnum + 1;
            coord(nnum, 1) = coordx;
            coord(nnum, 2) = coordy;
            coord(nnum, 3) = coordz;
        end
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
        idist = sqrt((coord(j, 1) - coord(i, 1))^2+(coord(j, 2) - coord(i, 2))^2+(coord(j, 3) - coord(i, 3))^2);
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
bc = 12.0 * emod / (pi * (delta^4));
% bc: Bond constant

disp = zeros(totnode, 3);
% disp: displacement of a material point, 1:x-coord, 2:y-coord, 3:z-coord
stendens = zeros(totnode, 3);
% stendens: strain energy of a material point, 1:loading 1, 2:loading 2, 3:loading 3
fncst = ones(totnode, 3);
% fncst: surface correction factors of a material point, 1:loading 1, 2:loading 2, 3:loading 3

% Loading 1
sedload1 = 0.6 * emod * 1.0e-6;
% sedload1: Strain energy density of a material point for the first loading condition

for i = 1:totnode
    disp(i, 1) = 0.001 * coord(i, 1);
    disp(i, 2) = 0.0;
    disp(i, 3) = 0.0;
end

for i = 1:totnode
    stendens(i, 1) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        %         idist = norm(coord(cnode, :)-coord(i, :));
        %         nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2+(coord(cnode, 3) - coord(i, 3))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2+(coord(cnode, 3) + disp(cnode, 3) - coord(i, 3) - disp(i, 3))^2);

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
sedload2 = 0.6 * emod * 1.0e-6;
% sedload2: Strain energy density of a material point for the second loading condition

for i = 1:totnode
    disp(i, 1) = 0.0;
    disp(i, 2) = 0.001 * coord(i, 2);
    disp(i, 3) = 0.0;
end

for i = 1:totnode
    stendens(i, 2) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        %         idist = norm(coord(cnode, :)-coord(i, :));
        %         nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2+(coord(cnode, 3) - coord(i, 3))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2+(coord(cnode, 3) + disp(cnode, 3) - coord(i, 3) - disp(i, 3))^2);

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

% Loading 3
sedload3 = 0.6 * emod * 1.0e-6;
% sedload3: Strain energy density of a material point for the third loading condition

for i = 1:totnode
    disp(i, 1) = 0.0;
    disp(i, 2) = 0.0;
    disp(i, 3) = 0.001 * coord(i, 3);
end

for i = 1:totnode
    stendens(i, 3) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        %         idist = norm(coord(cnode, :)-coord(i, :));
        %         nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2+(coord(cnode, 3) - coord(i, 3))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2+(coord(cnode, 3) + disp(cnode, 3) - coord(i, 3) - disp(i, 3))^2);

        if (idist <= delta - radij)
            fac = 1.0;
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij);
        else
            fac = 0.0;
        end

        stendens(i, 3) = stendens(i, 3) + 0.5 * 0.5 * bc * ((nlength - idist) / idist)^2 * idist * vol * fac;
    end
    % Calculation of surface correction factor in x direction
    % by finding the ratio of the analytical strain energy density value
    % to the strain energy density value obtained from PD Theory
    fncst(i, 3) = sedload3 / stendens(i, 3);
end

%% Initialization of displacements and velocities
vel = zeros(totnode, 3);
disp = zeros(totnode, 3);

%% Stable mass vector computation
dt = 1.0;
% dt: Time interval

massvec = zeros(totnode, 3);
% massvec: massvector for adaptive dynamic relaxation

for i = 1:totnode
    % 5 is a safety factor
    massvec(i, 1) = 0.25 * dt * dt * ((4.0 / 3.0) * pi * (delta)^3) * bc / dx;
    massvec(i, 2) = 0.25 * dt * dt * ((4.0 / 3.0) * pi * (delta)^3) * bc / dx;
    massvec(i, 3) = 0.25 * dt * dt * ((4.0 / 3.0) * pi * (delta)^3) * bc / dx;
end

%% Applied loading -  Right
appres = 200.0e6;
% appres: Applied pressure

bforce = zeros(totnode, 3);
% bforce: body load acting on a material point

% right
for i = 1:totint
    if (alflag(i, 1) == 1)
        bforce(i, 1) = appres / (dx);
    end
end

%% Time integration
nt = 1000;
% nt: Total number of time step
alpha = 23.0e-6;
% alpha: Coefficient of thermal expansion
dtemp = 0.0;
% dtemp: Temperature change
pratio = 1.0 / 4.0;
% pratio: Poisson's ratio

pforce = zeros(totnode, 3);
% pforce: total peridynamic force acting on a material point
pforceold = zeros(totnode, 3);
% pforceold: total peridynamic force acting on a material point in the previous time step 1:x-coord, 2:y-coord
acc = zeros(totnode, 3);
% acc: acceleration of a material point

velhalf = zeros(totnode, 3);
velhalfold = zeros(totnode, 3);
% vel: velocity of a material point

coord_disp_pd_ntbt = zeros(totint, 6);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
horiCnt = 0;
horizontal_dispsbt = zeros(ndivx, 9);
% Peridynamic displacement and Analytical displacement of points at y = 0;
vertiCnt = 0;
vertical_dispsbt = zeros(ndivy, 9);
% Peridynamic displacement and Analytical displacement of points at x = 0;
transCnt = 0;
transverse_dispsbt = zeros(ndivz, 9);

steady_check = zeros(nt, 4);


cn = 0.0;
cn1 = 0.0;
cn2 = 0.0;
for tt = 1:nt
    fprintf("%d/%d\n", tt, nt);

    for i = 1:totint
        pforce(i, 1) = 0.0;
        pforce(i, 2) = 0.0;
        pforce(i, 3) = 0.0;
        for j = 1:numfam(i, 1)
            cnode = nodefam(pointfam(i, 1)+j-1, 1);
            %             idist = norm(coord(cnode, :)-coord(i, :));
            %             nlength = norm((coord(cnode, :) + disp(cnode, :))-(coord(i, :) + disp(i, :)));
            idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2+(coord(cnode, 3) - coord(i, 3))^2);
            nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2+(coord(cnode, 3) + disp(cnode, 3) - coord(i, 3) - disp(i, 3))^2);

            % Volume correction
            if (idist <= delta - radij)
                fac = 1.0;
            elseif (idist <= delta + radij)
                fac = (delta + radij - idist) / (2.0 * radij);
            else
                fac = 0.0;
            end

            % Determination of the surface correction between two material points
            if (abs(coord(cnode, 3)-coord(i, 3)) <= 1.0e-10)
                if (abs(coord(cnode, 2)-coord(i, 2)) <= 1.0e-10)
                    theta = 0.0;
                elseif (abs(coord(cnode, 1)-coord(i, 1)) <= 1.0e-10)
                    theta = 90.0 * pi / 180.0;
                else
                    theta = atan(abs(coord(cnode, 2)-coord(i, 2))/abs(coord(cnode, 1)-coord(i, 1)));
                end
                phi = 90.0 * pi / 180.0;

                scx = (fncst(i, 1) + fncst(cnode, 1)) / 2.0;
                scy = (fncst(i, 2) + fncst(cnode, 2)) / 2.0;
                scz = (fncst(i, 3) + fncst(cnode, 3)) / 2.0;
                scr = 1.0 / (((cos(theta) * sin(phi))^2 / (scx)^2) + ((sin(theta) * sin(phi))^2 / (scy)^2) + ((cos(phi))^2 / (scz)^2));
                scr = sqrt(scr);
            elseif ((abs(coord(cnode, 1)-coord(i, 1)) <= 1.0e-10) && (abs(coord(cnode, 2)-coord(i, 2)) <= 1.0e-10))
                scz = (fncst(i, 3) + fncst(cnode, 3)) / 2.0;
                scr = scz;
            else
                theta = atan(abs(coord(cnode, 2)-coord(i, 2))/abs(coord(cnode, 1)-coord(i, 1)));
                phi = acos(abs(coord(cnode, 3)-coord(i, 3))/idist);

                scx = (fncst(i, 1) + fncst(cnode, 1)) / 2.0;
                scy = (fncst(i, 2) + fncst(cnode, 2)) / 2.0;
                scz = (fncst(i, 3) + fncst(cnode, 3)) / 2.0;
                scr = 1.0 / (((cos(theta) * sin(phi))^2 / (scx)^2) + ((sin(theta) * sin(phi))^2 / (scy)^2) + ((cos(phi))^2 / (scz)^2));
                scr = sqrt(scr);
            end

            % Calculation of the peridynamic force in x direction acting on a material point i due to a material point j
            dforce1 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1)) / nlength;
            dforce2 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2)) / nlength;
            dforce3 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (coord(cnode, 3) + disp(cnode, 3) - coord(i, 3) - disp(i, 3)) / nlength;

            pforce(i, 1) = pforce(i, 1) + dforce1;
            pforce(i, 2) = pforce(i, 2) + dforce2;
            pforce(i, 3) = pforce(i, 3) + dforce3;
        end
    end

    % Adaptive dynamic relaxation ⬇⬇⬇

    for i = 1:totint
        if (velhalfold(i, 1) ~= 0.0)
            cn1 = cn1 - disp(i, 1) * disp(i, 1) * (pforce(i, 1) / massvec(i, 1) - pforceold(i, 1) / massvec(i, 1)) / (dt * velhalfold(i, 1));
        end
        if (velhalfold(i, 2) ~= 0.0)
            cn1 = cn1 - disp(i, 2) * disp(i, 2) * (pforce(i, 2) / massvec(i, 2) - pforceold(i, 2) / massvec(i, 2)) / (dt * velhalfold(i, 2));
        end
        if (velhalfold(i, 3) ~= 0.0)
            cn1 = cn1 - disp(i, 3) * disp(i, 3) * (pforce(i, 3) / massvec(i, 3) - pforceold(i, 3) / massvec(i, 3)) / (dt * velhalfold(i, 3));
        end
        cn2 = cn2 + disp(i, 1) * disp(i, 1);
        cn2 = cn2 + disp(i, 2) * disp(i, 2);
        cn2 = cn2 + disp(i, 3) * disp(i, 3);
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

    for i = 1:totint
        % Integrate acceleration over time.
        if (tt == 1)
            velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
            velhalf(i, 2) = 1.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2)) / 2.0;
            velhalf(i, 3) = 1.0 * dt / massvec(i, 3) * (pforce(i, 3) + bforce(i, 3)) / 2.0;
        else
            velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
            velhalf(i, 2) = ((2.0 - cn * dt) * velhalfold(i, 2) + 2.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2))) / (2.0 + cn * dt);
            velhalf(i, 3) = ((2.0 - cn * dt) * velhalfold(i, 3) + 2.0 * dt / massvec(i, 3) * (pforce(i, 3) + bforce(i, 3))) / (2.0 + cn * dt);
        end

        vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
        vel(i, 2) = 0.5 * (velhalfold(i, 2) + velhalf(i, 2));
        vel(i, 3) = 0.5 * (velhalfold(i, 3) + velhalf(i, 3));

        disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;
        disp(i, 2) = disp(i, 2) + velhalf(i, 2) * dt;
        disp(i, 3) = disp(i, 3) + velhalf(i, 3) * dt;

        velhalfold(i, 1) = velhalf(i, 1);
        velhalfold(i, 2) = velhalf(i, 2);
        velhalfold(i, 3) = velhalf(i, 3);
        pforceold(i, 1) = pforce(i, 1);
        pforceold(i, 2) = pforce(i, 2);
        pforceold(i, 3) = pforce(i, 3);
    end

    % Adaptive dynamic relaxation ⬆⬆⬆

    if (tt == nt)
        for i = 1:totint
            coord_disp_pd_ntbt(i, 1:6) = [coord(i, 1), coord(i, 2), coord(i, 3), disp(i, 1), disp(i, 2), disp(i, 3)];
            if ((abs(coord(i, 2)-(dx / 2.0)) <= 1.0e-8) && (abs(coord(i, 3)-(dx / 2.0)) <= 1.0e-8))
                horiCnt = horiCnt + 1;
                horizontal_dispsbt(horiCnt, 1:9) = [coord(i, 1:3), ...
                    disp(i, 1:3), ...
                    0.001 * coord(i, 1), -1.0 * 0.001 * pratio * coord(i, 2), -1.0 * 0.001 * pratio * coord(i, 3)];
            end
            if ((abs(coord(i, 1)-(length / 2.0 + dx / 2.0)) <= 1.0e-8) && (abs(coord(i, 3)-(dx / 2.0)) <= 1.0e-8))
                vertiCnt = vertiCnt + 1;
                vertical_dispsbt(vertiCnt, 1:9) = [coord(i, 1:3), ...
                    disp(i, 1:3), ...
                    0.001 * coord(i, 1), -1.0 * 0.001 * pratio * coord(i, 2), -1.0 * 0.001 * pratio * coord(i, 3)];
            end
            if ((abs(coord(i, 1)-(length / 2.0 + dx / 2.0)) <= 1.0e-8) && (abs(coord(i, 2)-(dx / 2.0)) <= 1.0e-8))
                transCnt = transCnt + 1;
                transverse_dispsbt(transCnt, 1:9) = [coord(i, 1:3), ...
                    disp(i, 1:3), ...
                    0.001 * coord(i, 1), -1.0 * 0.001 * pratio * coord(i, 2), -1.0 * 0.001 * pratio * coord(i, 3)];
            end
        end
    end

    steady_check(tt, 1:4) = [tt, disp(7770, 1:3)];
end

%% plot
colormap jet;
factor = 1000;
subplot(241)
scatter3(coord_disp_pd_ntbt(:, 1)+factor*coord_disp_pd_ntbt(:, 4), ...
    coord_disp_pd_ntbt(:, 2)+factor*coord_disp_pd_ntbt(:, 5), ...
    coord_disp_pd_ntbt(:, 3)+factor*coord_disp_pd_ntbt(:, 6), [], ...
    sqrt(coord_disp_pd_ntbt(:, 4).^2+coord_disp_pd_ntbt(:, 5).^2+coord_disp_pd_ntbt(:, 6).^2), "filled")
subplot(242)
scatter3(coord_disp_pd_ntbt(:, 1)+factor*coord_disp_pd_ntbt(:, 4), ...
    coord_disp_pd_ntbt(:, 2)+factor*coord_disp_pd_ntbt(:, 5), ...
    coord_disp_pd_ntbt(:, 3)+factor*coord_disp_pd_ntbt(:, 6), [], ...
    coord_disp_pd_ntbt(:, 4), "filled")
subplot(245)
scatter3(coord_disp_pd_ntbt(:, 1)+factor*coord_disp_pd_ntbt(:, 4), ...
    coord_disp_pd_ntbt(:, 2)+factor*coord_disp_pd_ntbt(:, 5), ...
    coord_disp_pd_ntbt(:, 3)+factor*coord_disp_pd_ntbt(:, 6), [], ...
    coord_disp_pd_ntbt(:, 5), "filled")
subplot(246)
scatter3(coord_disp_pd_ntbt(:, 1)+factor*coord_disp_pd_ntbt(:, 4), ...
    coord_disp_pd_ntbt(:, 2)+factor*coord_disp_pd_ntbt(:, 5), ...
    coord_disp_pd_ntbt(:, 3)+factor*coord_disp_pd_ntbt(:, 6), [], ...
    coord_disp_pd_ntbt(:, 6), "filled")
subplot(243)
plot(steady_check(:, 1), steady_check(:, 2), steady_check(:, 1), steady_check(:, 3), steady_check(:, 1), steady_check(:, 4))
subplot(244)
plot(horizontal_dispsbt(:, 1), horizontal_dispsbt(:, 4), horizontal_dispsbt(:, 1), horizontal_dispsbt(:, 7))
subplot(247)
plot(vertical_dispsbt(:, 2), vertical_dispsbt(:, 5), vertical_dispsbt(:, 2), vertical_dispsbt(:, 8))
subplot(248)
plot(transverse_dispsbt(:, 3), transverse_dispsbt(:, 6), transverse_dispsbt(:, 3), transverse_dispsbt(:, 9))
