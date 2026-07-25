# Condition

A = 10;                    # radius
T = 60;                    # period
E = 0.00005;               # external force

sr = 100;                  # IMU sampling rate
dt = 1 / sr;               # simulation step

a = 0.1 * 4 * pi / T^2;    # angular acceleration
z = [A, 0, 0, 0, 0];       # initial state

t = 0:dt:(sqrt(10) * T);   # time
s = size(t)(2)

sigat = 700e-6 * 9.81;     # tang. accel. noise density
sigan = 700e-6 * 9.81;     # norm. accel. noise density
sigw = deg2rad(1);         # gyro noise density

sigats = sigat / sqrt(dt); # tang. accel. sampling standard deviation
sigans = sigan / sqrt(dt); # norm. accel. sampling standard deviation
sigws = sigw / sqrt(dt);   # gyro sampling standard deviation

# Reference trajectory

rw = a .* t;
ra = a .* t.^2 ./ 2;
rx = zeros(1, s);
ry = zeros(1, s);

rat = (a * A) .* ones(1, s);
ran = rw.^2 .* A;
rax = -rat .* sin(ra) - ran .* cos(ra);
ray = rat .* cos(ra) - ran .* sin(ra);

rx(1, 1) = z(1, 1);
ry(1, 1) = z(1, 2);
rvx = z(1, 3);
rvy = z(1, 4);

for k = 2:s
  rvx = rvx + rax(1, k) * dt;
  rvy = rvy + ray(1, k) * dt;
  rx(1, k) = rx(1, k - 1) + rvx * dt;
  ry(1, k) = ry(1, k - 1) + rvy * dt;
end

# Actual trajectory

yd = -E .* (t.^2 / 2);     # Y-disturbance due external force

aw = rw;
aa = ra;
ax = rx;
ay = ry + yd;

aat = rat + -E .* cos(aa);
aan = ran + E .* sin(aa);
aax = -aat .* sin(aa) - aan .* cos(aa);
aay = aat .* cos(aa) - aan .* sin(aa);

# DME

dmes = 5;           # DME sampling interval
dmei = dmes / dt;   # DME outage iterations
dmea = 0.1;         # DME accuracy (10cm)
dmex = 0;           # DME x-position
dmey = 0;           # DME y-position
dmed = zeros(1, s); # DME actual distances
coma = deg2rad(5);

for i = 1:s
  dmed(i) = sqrt((dmex - ax(i))^2 + (dmey - ay(i))^2);
end

dmem = dmed + dmea * randn(1, s);   # DME distance measurements
comm = aa + coma * randn(1, s);
R = diag([dmea^2]);                 # DME covariance
Rcom = diag([coma^2]);

# INS

iw = aw + sigws * randn(1, s);
iat = aat + sigats * randn(1, s);
ian = aan + sigans * randn(1, s);

# State setup

ppa = 0.00001;
paa = 0.00001;

# Pure
ix = zeros(1, s);
iy = zeros(1, s);
ia = zeros(1, s);

is = z;
ix(1, 1) = is(1, 1);
iy(1, 1) = is(1, 2);
ia(1, 1) = is(1, 5);

# Linearized Kalman
lx = zeros(1, s);
ly = zeros(1, s);
la = zeros(1, s);

ls = z;
lx(1, 1) = ls(1, 1);
ly(1, 1) = ls(1, 2);
la(1, 1) = ls(1, 5);

lQc = diag([sigat^2, sigan^2, sigw^2]);                        # continuous covariance
lesp = [0; 0; 0; 0; 0];                                        # prior error state
lPkp = diag([ppa^2, ppa^2, ppa^2, ppa^2, deg2rad(paa)^2]);     # prior error covariance

# EKF
ex = zeros(1, s);
ey = zeros(1, s);
ea = zeros(1, s);

es = z;
ex(1, 1) = es(1, 1);
ey(1, 1) = es(1, 2);
ea(1, 1) = es(1, 5);

eQc = diag([sigat^2, sigan^2, sigw^2]);                        # continuous covariance
eesp = [0; 0; 0; 0; 0];                                        # prior error state
ePkp = diag([ppa^2, ppa^2, ppa^2, ppa^2, deg2rad(paa)^2]);     # prior error covariance

# Runtime

J = [
  0, -1;
  1,  0
];

for k = 2:s
  # Pure

  # INS Integration

  is(1, 5) = is(1, 5) + iw(k) * dt;    # angle
  ia(k) = is(1, 5);

  iRot = [
    -sin(is(1, 5)), -cos(is(1, 5));
     cos(is(1, 5)), -sin(is(1, 5))
  ];

  iabf = [iat(k); ian(k)];             # body-frame acceleration
  ianf = iRot * iabf;                  # navigation-frame acceleration

  is(1, 3) = is(1, 3) + ianf(1) * dt;  # x-velocity
  is(1, 4) = is(1, 4) + ianf(2) * dt;  # y-velocity
  is(1, 1) = is(1, 1) + is(1, 3) * dt; # x
  is(1, 2) = is(1, 2) + is(1, 4) * dt; # y
  ix(1, k) = is(1, 1);
  iy(1, k) = is(1, 2);

  # Linearized Kalman

  # INS Integration

  ls(1, 5) = ls(1, 5) + iw(k) * dt;    # angle
  la(k) = ls(1, 5);

  lRot = [
    -sin(ra(1, k)), -cos(ra(1, k));
     cos(ra(1, k)), -sin(ra(1, k))
  ];

  labf = [
   iat(k);
   ian(k)
  ];                                   # body-frame acceleration
  lanf = lRot * labf;                  # navigation-frame acceleration

  ls(1, 3) = ls(1, 3) + lanf(1) * dt;  # x-velocity
  ls(1, 4) = ls(1, 4) + lanf(2) * dt;  # y-velocity
  ls(1, 1) = ls(1, 1) + ls(1, 3) * dt; # x
  ls(1, 2) = ls(1, 2) + ls(1, 4) * dt; # y
  lx(1, k) = ls(1, 1);
  ly(1, k) = ls(1, 2);

  lrefx = rx(1, k);
  lrefy = ry(1, k);
  lrefa = ra(1, k);

  # Filter

  ld0 = (lrefx - dmex) / sqrt((dmex - lrefx)^2 + (dmey - lrefy)^2);
  ld1 = (lrefy - dmey) / sqrt((dmex - lrefx)^2 + (dmey - lrefy)^2);

  lHk = [ld0, ld1, 0, 0, 0];
  lHcomk = [0, 0, 0, 0, 1];

  lgam = lRot * J * labf;
  lbet = -lRot;

  lFk = [
    zeros(2, 2), eye(2), zeros(2, 1);
    zeros(2, 4), lgam;
    zeros(1, 5)
  ];

  lGk = [
    zeros(2, 3);
    lbet, zeros(2, 1);
    0, 0, -1
  ];

  lCA = [
    -lFk, lGk * lQc * lGk';
    zeros(5, 5), lFk'
  ];

  lCB = expm(lCA);

  lPhik = lCB(6:10,6:10)';
  lQk = lPhik * lCB(1:5, 6:10);

  if (mod(k, dmei) == 0)
    # DME incorporation

    lKk = lPkp * lHk' * inv(lHk * lPkp * lHk' + R);
    lPk = (eye(5) - lKk * lHk) * lPkp;
    lzk = dmem(k) - sqrt((dmex - lrefx)^2 + (dmey - lrefy)^2);
    les = lesp + lKk * (lzk - lHk * lesp);

    ls = ls + les';
    les = zeros(5, 1);

    lesp = lPhik * les;
    lPkp = lPhik * lPk * lPhik' + lQk;
  else
    # DME outage

    lesp = lPhik * lesp;
    lPkp = lPhik * lPkp * lPhik' + lQk;
  endif

  # EKF

  # INS Integration

  es(1, 5) = es(1, 5) + iw(k) * dt;    # angle
  ea(k) = es(1, 5);

  eRot = [
    -sin(es(1, 5)), -cos(es(1, 5));
     cos(es(1, 5)), -sin(es(1, 5))
  ];

  eabf = [
   iat(k);
   ian(k)
  ];                                   # body-frame acceleration
  eanf = eRot * eabf;                  # navigation-frame acceleration

  es(1, 3) = es(1, 3) + eanf(1) * dt;  # x-velocity
  es(1, 4) = es(1, 4) + eanf(2) * dt;  # y-velocity
  es(1, 1) = es(1, 1) + es(1, 3) * dt; # x
  es(1, 2) = es(1, 2) + es(1, 4) * dt; # y
  ex(1, k) = es(1, 1);
  ey(1, k) = es(1, 2);

  erefx = es(1, 1);
  erefy = es(1, 2);
  erefa = es(1, 5);

  # Filter

  ed0 = (erefx - dmex) / sqrt((dmex - erefx)^2 + (dmey - erefy)^2);
  ed1 = (erefy - dmey) / sqrt((dmex - erefx)^2 + (dmey - erefy)^2);

  eHk = [ed0, ed1, 0, 0, 0];
  eHcomk = [0, 0, 0, 0, 1];

  egam = eRot * J * eabf;
  ebet = -eRot;

  eFk = [
    zeros(2, 2), eye(2), zeros(2, 1);
    zeros(2, 4), egam;
    zeros(1, 5)
  ];

  eGk = [
    zeros(2, 3);
    ebet, zeros(2, 1);
    0, 0, -1
  ];

  eCA = [
    -eFk, eGk * eQc * eGk';
    zeros(5, 5), eFk'
  ];

  eCB = expm(eCA);

  ePhik = eCB(6:10,6:10)';
  eQk = ePhik * eCB(1:5, 6:10);

  if (mod(k, dmei) == 0)
    # DME and Compass measurements available
    # Propagation
    eesp = ePhik * eesp;
    ePkp = ePhik * ePkp * ePhik' + eQk;

    # DME incorporation

    eKk = ePkp * eHk' * inv(eHk * ePkp * eHk' + R);
    ePk = (eye(5) - eKk * eHk) * ePkp;
    ezk = dmem(k) - sqrt((dmex - erefx)^2 + (dmey - erefy)^2);
    ees = eesp + eKk * (ezk - eHk * eesp);

    es = es + ees';       # update trajectory
    ees = zeros(5, 1);    # reset error state

    eesp = ees;           # update error state prior
    ePkp = ePk;           # update Pk prior

    # Compass incorporation

    eKk = ePkp * eHcomk' * inv(eHcomk * ePkp * eHcomk' + Rcom);
    ePk = (eye(5) - eKk * eHcomk) * ePkp;
    ezk = comm(k) - es(1, 5);
    ees = eesp + eKk * (ezk - eHcomk * eesp);

    es = es + ees';       # update trajectory
    ees = zeros(5, 1);    # reset error state

    eesp = ees;           # update error state prior
    ePkp = ePk;           # update Pk prior
  else
    # No DME and Compass measurements available
    # Propagation
    eesp = ePhik * eesp;
    ePkp = ePhik * ePkp * ePhik' + eQk;
  endif
end

figure(1);
plot(ax, ay, 'r-', ix, iy, 'b-', ex, ey, 'g-', lx, ly, 'y-');

figure(2);
plot(t, ax, 'r-', t, ix, 'b-');

figure(3);
plot(t, ay, 'r-', t, iy, 'b-');

figure(4);
plot(t, aa, 'r-', t, ia, 'b-');

figure(5);
plot(t, ax, 'r-', t, ex, 'g-');

figure(6);
plot(t, ay, 'r-', t, ey, 'g-');

figure(7);
plot(t, aa, 'r-', t, ea, 'g-');
