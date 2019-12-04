
%% Open/create new file 
% (Open or create new file for reading and writing in text mode. Discard any existing content)

fid = fopen(strcat(myFolder, filePrefix,'_simulationParameters.txt'), 'wt+');
fprintf(fid, 'date: %s\n', date);
fprintf(fid, '%s\n', '*** PARAMETERS ***');
fprintf(fid, 'flagDmLayout: %s\n', flagDmLayout);
fprintf(fid, 'ngsPhotometry: %s\n', char(ngsPhotometry));
fprintf(fid, 'sciPhotometry: %s\n', char(sciPhotometry));
fprintf(fid, 'nThresholded: %f\n', nThresholded);
fprintf(fid, 'controlBasis: %s\n', controlBasis);
fprintf(fid, 'pyrModulation: %f\n', pyrModulation);
fprintf(fid, 'startDelay: %f\n', startDelay);
fprintf(fid, 'exposureTime: %f\n', exposureTime);
fprintf(fid, 'gsMagnitude: %f\n', gsMagnitude);
fprintf(fid, 'readOutNoise: %f\n', readOutNoise);
fprintf(fid, 'photonNoise: %f\n', photonNoise);


fprintf(fid, '%s\n', '*** ATMOSPHERE ***');
fprintf(fid, 'fractionalR0: [%s]\n', strjoin(cellstr(num2str(fractionalR0(:))),', '));
fprintf(fid, 'altitude: [%s]\n', strjoin(cellstr(num2str(altitude(:))),', '));
fprintf(fid, 'windSpeed: [%s]\n', strjoin(cellstr(num2str(windSpeed(:))),', '));
fprintf(fid, 'windDirection: [%s]\n', strjoin(cellstr(num2str(windDirection(:))),', '));
fprintf(fid, 'L0: %f\n', L0);
fprintf(fid, 'r0: %f\n', r0);
fprintf(fid, 'seeing: %f\n', 0.98*0.5e-6/(constants.arcsec2radian*r0));
fprintf(fid, 'r0 Defined at: %f\n', atm.wavelength);

fprintf(fid, '%s\n', '*** TELESCOPE ***'); 
fprintf(fid, 'Number of lenslets (nL): %f\n', nL);
fprintf(fid, 'Number of pixels per lenslet (Npx): %f\n', nPx);
fprintf(fid, 'Resolution (pixels) on the pupil plane (nRes): %f\n', nRes);
fprintf(fid, 'Telescope primary mirror diameter (D): %f\n', D);
fprintf(fid, 'Lenslet pitch (d): %f\n', d);
fprintf(fid, 'WFS sampling time (samplingFreq): %f\n', samplingFreq);
fprintf(fid, 'Central obscuration ratio (obstructionRatio): %f\n', obstructionRatio);
fprintf(fid, 'Field of view (fieldOfViewInArcsec): %f\n', fieldOfViewInArcsec);

fprintf(fid, '%s\n', '*** KARHUNEN LOEVE FILE ***'); 
fprintf(fid, 'KLfilename: %s\n', KLfilename);

fprintf(fid, '%s\n', '*** WAVEFRONT SENSOR ***'); 

fprintf(fid, '%s\n', '*** DM ***'); 
% fprintf(fid, 'Fhr filename: %s\n', filename); 
% fprintf(fid, 'iFhr filename: %s\n', filenamei);


fprintf(fid, '%s\n', '*** CONTROLLER ***'); 
fprintf(fid, 'nDelay: %f\n', nDelay);

% integrator


fclose(fid);