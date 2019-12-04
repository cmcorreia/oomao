function quadDiff = fitPowerLaw(telemPSD, f, params)
 quadDiff = sum(abs(telemPSD - params(2)*f.^params(1)).^2);
 
