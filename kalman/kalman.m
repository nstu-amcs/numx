graphics_toolkit("gnuplot")

# Condition

A = 10;                    # radius
T = 60;                    # period
E = 0.0001;                # external force

sr = 100;                  # IMU sampling rate
dt = 1 / sr;               # simulation step

a = 0.1 * 4 * pi / T^2;    # angular acceleration
z = [A, 0, 0, 0, 0];       # initial state

t = 0:dt:(sqrt(10) * T);   # time
s = size(t)(2)

sigat = 400e-6 * 9.81;     # tang. accel. noise density
sigan = 400e-6 * 9.81;     # norm. accel. noise density
sigw = deg2rad(1);         # gyro noise density

sigats = sigat / sqrt(dt); # tang. accel. sampling standard deviation
sigans = sigan / sqrt(dt); # norm. accel. sampling standard deviation
sigws = sigw / sqrt(dt);   # gyro sampling standard deviation

# Reference trajectory

rw = a .* t;
ra = a .* t.^2 ./ 2;
rx = zeros(1, s);
ry = zeros(1, s);
rvx = zeros(1, s);
rvy = zeros(1, s);

rat = (a * A) .* ones(1, s);
ran = rw.^2 .* A;
rax = -rat .* sin(ra) - ran .* cos(ra);
ray = rat .* cos(ra) - ran .* sin(ra);

rx(1, 1) = z(1, 1);
ry(1, 1) = z(1, 2);
rvx(1, 1) = z(1, 3);
rvy(1, 1) = z(1, 4);

for k = 2:s
  rvx(1, k) = rvx(1, k - 1) + rax(1, k) * dt;
  rvy(1, k) = rvy(1, k - 1) + ray(1, k) * dt;
  rx(1, k) = rx(1, k - 1) + rvx(1, k) * dt;
  ry(1, k) = ry(1, k - 1) + rvy(1, k) * dt;
end

# Actual trajectory

yd = -E .* (t.^2 / 2);     # Y-disturbance due external force

aw = rw;
aa = ra;
ax = rx;
ay = ry; # + yd;
avx = rvx;
avy = rvy;

aat = rat; # -E .* cos(aa);
aan = ran; # + E .* sin(aa);
aax = -aat .* sin(aa) - aan .* cos(aa);
aay = aat .* cos(aa) - aan .* sin(aa);

# DME

dmes = 5;             # DME sampling interval
dmei = dmes / dt;     # DME outage iterations
dmea = 0.1;           # DME accuracy (10cm)

dme1x = 10;           # DME1 x-position
dme1y = 10;           # DME1 y-position
dme1d = zeros(1, s);  # DME1 actual distances

dme2x = -10;          # DME1 x-position
dme2y = -10;          # DME1 y-position
dme2d = zeros(1, s);  # DME1 actual distances

for i = 1:s
  dme1d(i) = sqrt((dme1x - ax(i))^2 + (dme1y - ay(i))^2);
  dme2d(i) = sqrt((dme2x - ax(i))^2 + (dme2y - ay(i))^2);
end

dme1m = dme1d + dmea * randn(1, s);   # DME1 distance measurements
dme2m = dme2d + dmea * randn(1, s);   # DME1 distance measurements

R = diag([dmea^2, dmea^2]);           # DME covariance
# R = diag([dmea^2]);

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
ivx = zeros(1, s);
ivy = zeros(1, s);
ia = zeros(1, s);

is = z;
ix(1, 1) = is(1, 1);
iy(1, 1) = is(1, 2);
ivx(1, 1) = is(1, 3);
ivy(1, 1) = is(1, 4);
ia(1, 1) = is(1, 5);

# Linearized Kalman
lx = zeros(1, s);
ly = zeros(1, s);
lvx = zeros(1, s);
lvy = zeros(1, s);
la = zeros(1, s);

ls = z;
lx(1, 1) = ls(1, 1);
ly(1, 1) = ls(1, 2);
lvx(1, 1) = ls(1, 3);
lvy(1, 1) = ls(1, 4);
la(1, 1) = ls(1, 5);

lQc = diag([sigat^2, sigan^2, sigw^2]);                        # continuous covariance
lesp = [0; 0; 0; 0; 0];                                        # prior error state
lPkp = diag([ppa^2, ppa^2, ppa^2, ppa^2, deg2rad(paa)^2]);     # prior error covariance

# EKF
ermsx = zeros(1, s);
ermsy = zeros(1, s);

ex = zeros(1, s);
ey = zeros(1, s);
evx = zeros(1, s);
evy = zeros(1, s);
ea = zeros(1, s);

es = z;
ex(1, 1) = es(1, 1);
ey(1, 1) = es(1, 2);
evx(1, 1) = es(1, 3);
evy(1, 1) = es(1, 4);
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
  ivx(1, k) = is(1, 3);
  ivy(1, k) = is(1, 4);

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
  lvx(1, k) = ls(1, 3);
  lvy(1, k) = ls(1, 4);

  lrefx = rx(1, k);
  lrefy = ry(1, k);

  # Filter

  ld10 = (lrefx - dme1x) / sqrt((dme1x - lrefx)^2 + (dme1y - lrefy)^2);
  ld11 = (lrefy - dme1y) / sqrt((dme1x - lrefx)^2 + (dme1y - lrefy)^2);
  ld20 = (lrefx - dme2x) / sqrt((dme2x - lrefx)^2 + (dme2y - lrefy)^2);
  ld21 = (lrefy - dme2y) / sqrt((dme2x - lrefx)^2 + (dme2y - lrefy)^2);

  lHk = [
    ld10, ld11, 0, 0, 0;
    ld20, ld21, 0, 0, 0
  ];

  # lHk = [ld10, ld11, 0, 0, 0];

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
  ] * dt;

  lCB = expm(lCA);

  lPhik = lCB(6:10,6:10)';
  lQk = lPhik * lCB(1:5, 6:10);

  if (mod(k, dmei) == 0)
    # DME incorporation

    lKk = lPkp * lHk' * inv(lHk * lPkp * lHk' + R);

    lzk = [
      dme1m(k) - sqrt((dme1x - lrefx)^2 + (dme1y - lrefy)^2);
      dme2m(k) - sqrt((dme2x - lrefx)^2 + (dme2y - lrefy)^2)
    ];

    # lzk = [dme1m(k) - sqrt((dme1x - lrefx)^2 + (dme1y - lrefy)^2)];

    les = lesp + lKk * (lzk - lHk * lesp);
    lPk = lPkp - lKk * (lHk * lPkp * lHk' + R) * lKk';

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
  evx(1, k) = es(1, 3);
  evy(1, k) = es(1, 4);

  erefx = es(1, 1);
  erefy = es(1, 2);

  # Filter

  ed10 = (erefx - dme1x) / sqrt((dme1x - erefx)^2 + (dme1y - erefy)^2);
  ed11 = (erefy - dme1y) / sqrt((dme1x - erefx)^2 + (dme1y - erefy)^2);
  ed20 = (erefx - dme2x) / sqrt((dme2x - erefx)^2 + (dme2y - erefy)^2);
  ed21 = (erefy - dme2y) / sqrt((dme2x - erefx)^2 + (dme2y - erefy)^2);

  eHk = [
    ed10, ed11, 0, 0, 0;
    ed20, ed21, 0, 0, 0
  ];

  # eHk = [ed10, ed11, 0, 0, 0];

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
  ] * dt;

  eCB = expm(eCA);

  ePhik = eCB(6:10,6:10)';
  eQk = ePhik * eCB(1:5, 6:10);

  if (mod(k, dmei) == 0)
    # DME available -> incorporation

    eKk = ePkp * eHk' * inv(eHk * ePkp * eHk' + R);            # gain

    ezk = [
      dme1m(k) - sqrt((dme1x - erefx)^2 + (dme1y - erefy)^2);
      dme2m(k) - sqrt((dme2x - erefx)^2 + (dme2y - erefy)^2)
    ];

    # ezk = [dme1m(k) - sqrt((dme1x - erefx)^2 + (dme1y - erefy)^2)];

    ees = eesp + eKk * (ezk - eHk * eesp);                     # update error-state estimate
    ePk = ePkp - eKk * (eHk * ePkp * eHk' + R) * eKk';         # updated estimate cov

    es = es + ees';    # INS state update
    ees = zeros(5, 1); # error-state reset

    # Projection
    eesp = ePhik * ees;                  # next prior error-state
    ePkp = ePhik * ePk * ePhik' + eQk;   # next prior error-state cov
  else
    # No DME available -> projection only

    ees = eesp;
    ePk = ePkp;

    # Projection
    eesp = ePhik * ees;                  # next prior error-state
    ePkp = ePhik * ePk * ePhik' + eQk;   # next prior error-state cov
  endif
end

grid on;

subplot(2, 3, 1);
plot(ax, ay, 'k--', ix, iy, 'k:', ex, ey, 'k-', lx, ly, 'k-');
title('Trajectory');

subplot(2, 3, 2);
plot(t, ax, 'k--', t, ix, 'k:', t, ex, 'k-', t, lx, 'k-');
title('Coordinate X');

subplot(2, 3, 3);
plot(t, avx, 'k--', t, ivx, 'k:', t, evx, 'k-', t, lvx, 'k-');
title('Velocity X');

subplot(2, 3, 4);
plot(t, aa, 'k--', t, ia, 'k:', t, ea, 'k-', t, la, 'k-');
title('Angle');

subplot(2, 3, 5);
plot(t, ay, 'k--', t, iy, 'k:', t, ey, 'k-', t, ly, 'k-');
title('Coordinate Y');

subplot(2, 3, 6);
plot(t, avy, 'k--', t, ivy, 'k:', t, evy, 'k-', t, lvy, 'k-');
title('Velocity Y');

difx = sum((ax - ex).^2);
dify = sum((ay - ey).^2);
dift = difx + dify;

gamx = sqrt(difx / s);
gamy = sqrt(dify / s);
gamt = sqrt(dift / (2 * s));
