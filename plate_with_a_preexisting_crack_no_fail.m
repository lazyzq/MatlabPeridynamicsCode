%%
clear;
clc;

%% Specification of the locations of material points
length = 0.05;
% length: Total length of the bar
width = 0.05;
% width: Total width of the plate

ndivx = 500;
% ndivx: Number of divisions in x direction - except boundary region
ndivy = 500;
% ndivy: Number of divisions in y direction - except boundary region

dx = length / ndivx;
% dx: Spacing between material points
nbnd = 3;
% nbnd: Number of divisions in the boundary region
totnode = ndivx * (ndivy + 2*nbnd);
% totnode: Total number of material points

coord = zeros(totnode, 2);
% coord: Material point locations

nnum = 0;
% nnum: Material point number

% Material points of the internal region
for i = 1:ndivy
    for j = 1:ndivx
        nnum = nnum + 1;
        coord(nnum, 1) = (-1.0*length/2.0) + (dx/2.0) + (j - 1)*dx;
        coord(nnum, 2) = (-1.0*width/2.0) + (dx/2.0) + (i - 1)*dx;
    end
end

totint = nnum;

% Material points of the boundary region - bottom
for i = 1:nbnd
    for j = 1: ndivx
        nnum = nnum + 1;
        coord(nnum, 1) = -1.0/2.0*length + (dx/2.0) + (j - 1)*dx;
        coord(nnum, 2) = -1.0/2.0*width - (dx/2.0) - (i - 1)*dx;
    end
end

totbottom = nnum;

% Material points of the boundary region - top
for i = 1: nbnd
    for j = 1: ndivx
        nnum = nnum + 1;
        coord(nnum, 1) = -1.0/2.0*length + (dx/2.0) + (j - 1)*dx;
        coord(nnum, 2) = 1.0/2.0*width + (dx/2.0) + (i - 1)*dx;
    end
end

tottop = nnum;
totnode = nnum;


%% Initialization of fail flag array
maxfam = 100;
% maxfam: Maximum number of material points inside a horizon of a material point

fail = ones(totnode, maxfam);

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

%% Definition of the crack surface
crlength = 0.01;
% crlength: Crack length

% PD bonds penetrating through the crack surface are broken
for i = 1: totnode
    for j = 1: numfam(i, 1)
        cnode = nodefam(pointfam(i, 1) + j - 1, 1);
        if ((coord(cnode, 2) > 0.0) && (coord(i, 2) < 0.0))
            if ((abs(coord(i, 1)) - (crlength/2.0)) <= 1.0e-10)
                fail(i, j) = 0;
            elseif ((abs(coord(cnode, 1)) - (crlength/2.0)) <= 1.0e-10)
                fail(i, j) = 0;
            end
        elseif ((coord(i, 2) > 0.0) && (coord(cnode, 2) < 0.0))
            if ((abs(coord(i, 1)) - (crlength/2.0)) <= 1.0e-10)
                fail(i, j) = 0;
            elseif ((abs(coord(cnode, 1)) - (crlength/2.0)) <= 1.0e-10)
                fail(i, j) = 0;
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

dens = 8000.0;
% dens: Density
emod = 192.0e9;
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
dt = 0.8*sqrt(2.0*dens*dx/(pi*delta^2*dx*bc));
% dt: Time interval

% massvec = zeros(totnode, 2);
% % massvec: massvector for adaptive dynamic relaxation
% 
% for i = 1:totnode
%     % 5 is a safety factor
%     massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc / dx;%  * 5.0;
%     massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc / dx;%  * 5.0;
% end

%% Applied loading - Left and Right
appres = 200.0e6;
% appres: Applied pressure

bforce = zeros(totnode, 2);
% bforce: body load acting on a material point


%% Time integration
nt = 1250;
% nt: Total number of time step
alpha = 23.0e-6;
% alpha: Coefficient of thermal expansion
dtemp = 0.0;
% dtemp: Temperature change
pratio = 1.0 / 3.0;
% pratio: Poisson's ratio
scr0 = 1.0;
% scr0: Critical stretch

dmg = zeros(totint, 1);

pforce = zeros(totnode, 2);
% pforce: total peridynamic force acting on a material point
pforceold = zeros(totnode, 2);
% pforceold: total peridynamic force acting on a material point in the previous time step 1:x-coord, 2:y-coord
acc = zeros(totnode, 2);
% acc: acceleration of a material point

velhalf = zeros(totnode, 2);
velhalfold = zeros(totnode, 2);
% vel: velocity of a material point

coord_disp_pd_750_pwc = zeros(totint, 4);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
coord_disp_pd_1000_pwc = zeros(totint, 4);
coord_disp_pd_1250_pwc = zeros(totint, 4);

steady_check = zeros(nt, 3);

cn = 0.0;
cn1 = 0.0;
cn2 = 0.0;
for tt = 1:nt
    fprintf("%d/%d\n", tt, nt);
    ctime = tt * dt;

    % Application of boundary conditions at the top and bottom edges
    for i = (totint + 1): totbottom
        vel(i, 2) = -20.0;
        disp(i, 2) = -20.0*tt*dt;
    end

    for i = (totbottom + 1): tottop
        vel(i, 2) = 20.0;
        disp(i, 2) = 20.0*tt*dt;
    end


    for i = 1:totint
        dmgpar1 = 0.0;
        dmgpar2 = 0.0;
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

            if (fail(i, j) == 1)
                dforce1 = bc*(nlength - idist)/idist*vol*scr*fac*(coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))/nlength;
                dforce2 = bc*(nlength - idist)/idist*vol*scr*fac*(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))/nlength;
            else
                    dforce1 = 0.0;
                    dforce2 = 0.0;
            end
            pforce(i, 1) = pforce(i, 1) + dforce1;
            pforce(i, 2) = pforce(i, 2) + dforce2;

            % Definition of a no-fail zone
            if (abs((nlength - idist)/idist) > scr0)
                if (abs(coord(i, 2)) <= (length/4.0))
                    fail(i, j) = 0;
                end
            end

            dmgpar1 = dmgpar1 + fail(i, j)*vol*fac;
            dmgpar2 = dmgpar2 + vol*fac;
        end
        % Calculation of the damage parameter 
        dmg(i, 1) = 1.0 - dmgpar1/dmgpar2;
    end
    for i = 1: totint
        % Calculation of acceleration of material point i
        acc(i, 1) = (pforce(i, 1) + bforce(i, 1))/dens;
        acc(i, 2) = (pforce(i, 2) + bforce(i, 2))/dens;
        % Calculation of velocity of material point i
        % by integrating the acceleration of material point i
        vel(i, 1) = vel(i, 1) + acc(i, 1)*dt;
        vel(i, 2) = vel(i, 2) + acc(i, 2)*dt;
        % Calculation of displacement of material point i
        % by integrating the velocity of material point i
        disp(i, 1) = disp(i, 1) + vel(i, 1)*dt;
        disp(i, 2) = disp(i, 2) + vel(i, 2)*dt;
    end

    for i = (totint + 1): totbottom
        acc(i, 1) = (pforce(i, 1) + bforce(i, 1))/dens;
        vel(i, 1) = vel(i, 1) + acc(i, 1)*dt;
        disp(i, 1) = disp(i, 1) + vel(i, 1)*dt;
    end

    for i = (totbottom + 1): tottop
        acc(i, 1) = (pforce(i, 1) + bforce(i, 1))/dens;
        vel(i, 1) = vel(i, 1) + acc(i, 1)*dt;
        disp(i, 1) = disp(i, 1) + vel(i, 1)*dt;
    end

    if (tt == 750)
        for i = 1:totint
           coord_disp_pd_750_pwc(i, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
        end
    elseif (tt == 1000)
        for i = 1:totint
           coord_disp_pd_1000_pwc(i, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
        end
    elseif (tt == 1250)
        for i = 1:totint
           coord_disp_pd_1250_pwc(i, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
        end
    end

    steady_check(tt, 1:3) = [tt, disp(3788, 1), disp(3788, 2)];
end

%% plot
colormap jet;
subplot(221)
plot(steady_check(:, 1), steady_check(:, 2), steady_check(:, 1), steady_check(:, 3))
scale = 5;
subplot(222)
scatter(coord_disp_pd_750_pwc(:, 1)+scale*coord_disp_pd_750_pwc(:, 3), coord_disp_pd_750_pwc(:, 2)+scale*coord_disp_pd_750_pwc(:, 4), [], sqrt(coord_disp_pd_750_pwc(:, 3).^2+coord_disp_pd_750_pwc(:, 4).^2), "filled")
subplot(223)
scatter(coord_disp_pd_1000_pwc(:, 1)+scale*coord_disp_pd_1000_pwc(:, 3), coord_disp_pd_1000_pwc(:, 2)+scale*coord_disp_pd_1000_pwc(:, 4), [], sqrt(coord_disp_pd_1000_pwc(:, 3).^2+coord_disp_pd_1000_pwc(:, 4).^2), "filled")
subplot(224)
scatter(coord_disp_pd_1250_pwc(:, 1)+scale*coord_disp_pd_1250_pwc(:, 3), coord_disp_pd_1250_pwc(:, 2)+scale*coord_disp_pd_1250_pwc(:, 4), [], sqrt(coord_disp_pd_1250_pwc(:, 3).^2+coord_disp_pd_1250_pwc(:, 4).^2), "filled")
