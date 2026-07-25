graphics_toolkit("gnuplot");

# Condition

A = 10;                   # radius
T = 60;                   # period
E = 0.005;                 # external force

sr = 100;                 # IMU sampling rate
dt = 1 / sr;              # simulation step

a = 0.1 * 4 * pi / T^2;   # angular acceleration
z = [A, 0, 0, 0, 0];      # initial state

t = 0:dt:(sqrt(10)*T);           # time
s = size(t)(2)

sigat = 400e-6 * 9.81;   # tang. accel. noise density
sigan = 400e-6 * 9.81;   # norm. accel. noise density
sigw = deg2rad(0.001);    # gyro noise density

sigats = sigat / sqrt(dt); # tang. accel. sampling standard deviation
sigans = sigan / sqrt(dt); # norm. accel. sampling standard deviation
sigws = sigw / sqrt(dt);   # gyro sampling standard deviation

# Reference trajectory

rw = a .* t;
ra = a .* t.^2 ./ 2;
rx = A .* cos(ra);
ry = A .* sin(ra);

rat = (a * A) .* ones(1, s);
ran = rw.^2 .* A;
rax = -rat .* sin(ra) - ran .* cos(ra);
ray = rat .* cos(ra) - ran .* sin(ra);

# Actual trajectory

yd = -E .* (t.^2 / 2); # Y-disturbance due external force

aw = rw;
aa = ra;
ax = rx;
ay = ry; # + yd;

eat = -E .* cos(aa);
ean = E .* sin(aa);
aat = rat; # + eat;
aan = ran; # + ean;
aax = -aat .* sin(ra) - aan .* cos(ra);
aay = aat .* cos(ra) - aan .* sin(ra);

# DME

dmes = 5;          # DME sampling interval
dmei = dmes / dt;  # DME outage iterations
dmea = 0.1;        # DME accuracy (10cm)
dmex = 0;          # DME x-position
dmey = 0;          # DME y-position

dmed = zeros(1, s); # DME actual distances

for i = 1:s
  dmed(i) = sqrt((dmex - ax(i))^2 + (dmey - ay(i))^2);
end

dmem = dmed + dmea * randn(1, s);   # DME distance measurements
R = diag([dmea^2]);                 # DME covariance

# INS

ix = zeros(1, s);
iy = zeros(1, s);
ia = zeros(1, s);

iw = aw + sigws * randn(1, s);
iat = aat + sigats * randn(1, s);
ian = aan + sigans * randn(1, s);

# Initial state setup

is = z;
ix(1, 1) = is(1, 1);
iy(1, 1) = is(1, 2);
ia(1, 1) = is(1, 5);

# Runtime

Qc = diag([sigat^2, sigan^2, sigw^2]);                        # continuous covariance
esp = [0; 0; 0; 0; 0];                                        # prior error state
Pkp = diag([0.01^2, 0.01^2, 0.01^2, 0.01^2, deg2rad(0.1)^2]); # prior error covariance

for i = 2:s
  # INS integration

  is(1, 5) = is(1, 5) + iw(i) * dt;     # angle
  ca = is(1, 5);
  ia(i) = ca;

  # Linearized Kalman reference

  #Rot = [
  #  -sin(ra(i)), -cos(ra(i));
  #   cos(ra(i)), -sin(ra(i))
  #];

  # EKF reference

  Rot = [
    -sin(ca), -cos(ca);
     cos(ca), -sin(ca)
  ];

  ##

  J = [
    0, -1;
    1,  0
  ];

  iabf = [iat(i); ian(i)];               # INS body-frame acceleration
  ianf = Rot * iabf;                     # INS navigation-frame acceleration

  is(1, 3) = is(1, 3) + ianf(1) * dt;  # x-velocity
  is(1, 4) = is(1, 4) + ianf(2) * dt;  # y-velocity

  is(1, 1) = is(1, 1) + is(1, 3) * dt; # x
  is(1, 2) = is(1, 2) + is(1, 4) * dt; # y

  ix(1, i) = is(1, 1);
  iy(1, i) = is(1, 2);

  # Linearized Kalman reference

  #refx = rx(1, i);
  #refy = ry(1, i);

  # EKF reference

  refx = is(1, 1);
  refy = is(1, 2);

  ##

  # Filter

  d0 = (refx - dmex) / sqrt((dmex - refx)^2 + (dmey - refy)^2);
  d1 = (refy - dmey) / sqrt((dmex - refx)^2 + (dmey - refy)^2);

  Hk = [d0, d1, 0, 0, 0];

  gam = Rot * J * iabf;
  bet = -Rot;

  Fk = [
    zeros(2, 2), eye(2), zeros(2, 1);
    zeros(2, 4), gam;
    zeros(1, 5)
  ];

  Gk = [
    zeros(2, 3);
    bet, zeros(2, 1);
    0, 0, -1
  ];

  CA = [
    -Fk, Gk * Qc * Gk';
    zeros(5, 5), Fk'
  ];

  CB = expm(CA);

  Phik = CB(6:10,6:10)';
  Qk = Phik * CB(1:5, 6:10);

  esp = Phik * esp;
  Pkp = Phik * Pkp * Phik' + Qk;
end

figure(1);
plot(rx, ry, 'r-', ix, iy, 'b-');

figure(2);
plot(t, rx, 'r-', t, ax, 'g-', t, ix, 'b-');

figure(3);
plot(t, ry, 'r-', t, ay, 'g-', t, iy, 'b-');

figure(4);
plot(t, ra, 'r-', t, aa, 'g-', t, ia, 'b-');
