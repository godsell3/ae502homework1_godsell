clear;
clc;
format longg;

%{
1/1/2017 - Day 1
12/31/2017 - Day 365
8/1/2017 - Day 213
1/31/2019 - Day 761
1/1/2017 - Day 1
7/31/2020 - Day 1308
6/1/2019 - Day 882
1/31/2022 - Day 1857
%}

% Inputs provided in the problem
r1I = [3.515868886595499e-2, -3.162046390773074, 4.493983111703389]; %au
v1I = [-2.317577766980901e-3, 9.843360903693031e-3, -1.541856855538041e-2]; %au/day
r2I = [7.249472033259724, 14.61063037906177, 14.24274452216359]; %au
v2I = [-8.241709369476881e-3, -1.156219024581502e-2, -1.317135977481448e-2]; %au/day
rE  = [-1.796136509111975e-1, 9.667949206859814e-1, -3.668681017942158e-5]; %au
vE  = [-1.720038360888334e-2, -3.211186197806460e-3, 7.927736735960840e-7]; %au/day

% Canonical units calculation (want DU=AU, TU=?days)
DU = 149597870.; %km
mu_sun = 1.327e11; %km^3/s^2
P = 2*pi*sqrt(DU^3/mu_sun); %period of Earth about sun
TU = P/(2*pi); %s
TU_days = TU/(3600*24); %days

% Convert input vectors to canonical units (AU=DU)
v1I = v1I*TU_days;
v2I = v2I*TU_days;
vE = vE*TU_days;

% Define number of steps
steps = 2000;

% Earth trajectory
dt = 365.25/TU_days; %TU in a year
% [rf, vf] = findFinal_rv(rE, vE, dt, 1);
[E_traj, E_tvals] = univ_propagator(rE, vE, dt, steps, 1);
dt = 1308/TU_days; %TU in Borisov launch window
[E_traj2, E_tvals2] = univ_propagator(rE, vE, dt, steps, 1);

% 1I trajectory
dt = 760/TU_days; %TU from 1/1/2017 to 1/31/2019
% [rf, vf] = findFinal_rv(r1I, v1I, dt, 1);
[I1_traj, I1_tvals] = univ_propagator(r1I, v1I, dt, steps, 1);
[~, split_idx] = min(abs(I1_tvals.*TU_days-212));
I1_traj = I1_traj(split_idx:end,:);
I1_tvals = I1_tvals(split_idx:end);

% 2I trajectory
dt = 365.25/TU_days*4; %TU in 4 years
% [rf, vf] = findFinal_rv(r2I, v2I, dt, 1);
[I2_traj, I2_tvals] = univ_propagator(r2I, v2I, dt, steps, 1);
[~, split_idx] = min(abs(I2_tvals.*TU_days-881));
I2_traj = I2_traj(split_idx:end,:);
I2_tvals = I2_tvals(split_idx:end);

% Plot trajectories
figure(1)
hold on
grid on
axis equal
box on
plot3(E_traj(:,1), E_traj(:,2), E_traj(:,3)) %plot Earth
plot3(E_traj(end,1),E_traj(end,2),E_traj(end,3),'.','MarkerSize',10) %plot Earth endpoint
plot3(0.,0.,0.,'.','MarkerSize',20) %plot sun
plot3(I1_traj(:,1), I1_traj(:,2), I1_traj(:,3)) %plot 1I
plot3(I1_traj(end,1),I1_traj(end,2),I1_traj(end,3),'.','MarkerSize',10) %plot 1I endpoint
plot3(I2_traj(:,1), I2_traj(:,2), I2_traj(:,3)) %plot 2I
plot3(I2_traj(end,1),I2_traj(end,2),I2_traj(end,3),'.','MarkerSize',10) %plot 2I endpoint
title('Earth, 1I, and 2I trajectories about the Sun')
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
view(20,15)
legend('Earth trajectory','Earth (12/31/2017)','Sun','1I trajectory','1I (1/31/2019)','2I trajectory','2I (1/31/2022)')


% Create datetime objects
t0 = datetime(2017,1,1);
te1 = datetime(2017,12,31);
t1I_1 = datetime(2017,8,1);
t1I_2 = datetime(2019,1,31);
te2 = datetime(2020,7,31);
t2I_1 = datetime(2019,6,1);
t2I_2 = datetime(2022,1,31);

% Retrieve pork chop data for 1I
[deltav_1I,deltav_1Ifb] = porkchopSolver(E_traj,E_tvals,I1_traj,I1_tvals,'prograde',1);
deltav_1I(deltav_1I==-1) = 0; %eliminate negative values
deltav_1Ifb(deltav_1Ifb==-1) = 0;
deltav_1I = deltav_1I*DU/TU; %convert to km/s
deltav_1Ifb = deltav_1Ifb*DU/TU;
deltav_1I(deltav_1I>50) = 0; %eliminate values >50 km/s
deltav_1Ifb(deltav_1Ifb>20) = 0; %eliminate values >20 km/s
deltav_1I = deltav_1I'; %transpose
deltav_1Ifb = deltav_1Ifb';
deltav_1I(deltav_1I==0) = NaN; %delete irrelevant values
deltav_1Ifb(deltav_1Ifb==0) = NaN;

% Plot pork chop data for 1I (rendezvous)
figure(2)
dep_range = linspace(t0, te1, length(deltav_1I(1,:)));
arr_range = linspace(t1I_1, t1I_2, length(deltav_1I(:,1)));
porkchop1Ir = surf(dep_range, arr_range, deltav_1I);
porkchop1Ir.EdgeColor = 'none';
xlabel('Departure Date')
xlim([t0 te1])
ylabel('Arrival Date')
ylim([t1I_1 t1I_2])
title('Prograde Rendez-vous with 1I/`Oumuamua from Earth')
colormap(jet)
a = colorbar;
a.Label.String = '\DeltaV';
a.Label.Rotation = 0;
a.Label.Position(1) = 3;
view(0,90)

% Plot pork chop data for 1I (flyby)
figure(3)
dep_range = linspace(t0, te1, length(deltav_1Ifb(1,:)));
arr_range = linspace(t1I_1, t1I_2, length(deltav_1Ifb(:,1)));
porkchop1Ir = surf(dep_range, arr_range, deltav_1Ifb);
porkchop1Ir.EdgeColor = 'none';
xlabel('Departure Date')
xlim([t0 datetime(2017,11,30)])
ylabel('Arrival Date')
ylim([t1I_1 datetime(2018,7,31)])
title('Prograde Fly-By with 1I/`Oumuamua from Earth')
colormap(jet)
a = colorbar;
a.Label.String = '\DeltaV';
a.Label.Rotation = 0;
a.Label.Position(1) = 3;
view(0,90)

% Retrieve pork chop data for 2I
[deltav_2I,deltav_2Ifb] = porkchopSolver(E_traj2,E_tvals2,I2_traj,I2_tvals,'prograde',1);
deltav_2I(deltav_2I==-1) = 0; %eliminate negative values
deltav_2Ifb(deltav_2Ifb==-1) = 0;
deltav_2I = deltav_2I*DU/TU; %convert to km/s
deltav_2Ifb = deltav_2Ifb*DU/TU;
deltav_2I(deltav_2I>60) = 0; %eliminate values >60 km/s
deltav_2Ifb(deltav_2Ifb>20) = 0; %eliminate values >20 km/s
deltav_2I = deltav_2I'; %transpose
deltav_2Ifb = deltav_2Ifb';
deltav_2I(deltav_2I==0) = NaN; %delete irrelevant values
deltav_2Ifb(deltav_2Ifb==0) = NaN;

% Plot pork chop data for 2I (rendezvous)
figure(4)
dep_range = linspace(t0, te2, length(deltav_2I(1,:)));
arr_range = linspace(t2I_1, t2I_2, length(deltav_2I(:,1)));
porkchop2Ir = surf(dep_range, arr_range, deltav_2I);
porkchop2Ir.EdgeColor = 'none';
xlabel('Departure Date')
xlim([t0 te2])
ylabel('Arrival Date')
ylim([t2I_1 t2I_2])
title('Prograde Rendez-vous with 2I/Borisov from Earth')
colormap(jet)
b = colorbar;
b.Label.String = '\DeltaV';
b.Label.Rotation = 0;
b.Label.Position(1) = 3;
view(0,90)

% Plot pork chop data for 2I (flyby)
figure(5)
dep_range = linspace(t0, te2, length(deltav_2Ifb(1,:)));
arr_range = linspace(t2I_1, t2I_2, length(deltav_2Ifb(:,1)));
porkchop1Ir = surf(dep_range, arr_range, deltav_2Ifb);
porkchop1Ir.EdgeColor = 'none';
xlabel('Departure Date')
xlim([t0 datetime(2019,1,31)])
ylabel('Arrival Date')
ylim([t2I_1 datetime(2020,8,31)])
title('Prograde Fly-By with 2I/Borisov from Earth')
colormap(jet)
a = colorbar;
a.Label.String = '\DeltaV';
a.Label.Rotation = 0;
a.Label.Position(1) = 3;
view(0,90)

% Convert the initial state vectors into orbital elements
[a_1I, e_1I, i_1I, RAAN_1I, w_1I, f_1I] = orbitalElements(r1I, v1I, 1);
a_1I = a_1I*DU; %km
[a_2I, e_2I, i_2I, RAAN_2I, w_2I, f_2I] = orbitalElements(r2I, v2I, 1);
a_2I = a_2I*DU; %km
