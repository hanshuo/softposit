imagePts = [612, 486, 567, 476, 441, 117, 145, 206, 234, 329];
imagePts = reshape(imagePts,[],2)

worldPts = [-3.75, 7.50, -3, 3, 0, 0, 0, 0, -5.0, 5.0, 2.25, -2.25, 0.5, 2.75, -2.0, -2.0, -0.75, -0.75];
worldPts = reshape(worldPts,[],3)

beta0 = 2.0E-4;
noiseStd = 10.0;

initRot = [0, 1, 0, -1, 0, 0, 0, 0, 1];
initRot = reshape(initRot,[],3)

initTrans = [0.0; 0.0; 30.0];
focalLength = 982.1;
center = [376, 240];

[rot, trans] = SoftPOSIT(imagePts, worldPts, beta0, noiseStd, initRot, initTrans, focalLength, center)

