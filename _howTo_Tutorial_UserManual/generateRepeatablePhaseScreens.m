%% EXAMPLE MATLAB RANDOM GENERATION

% Create a stream (~=seed)
rngStream = RandStream('mt19937ar');
initState = rngStream.State;
rand(rngStream,1,5)
rand(rngStream,1,5)
rand(rngStream,1,5)

%reset :: the rand function will output the same numbers as it did at the
%beginning
rngStream.State = initState;
rand(rngStream,1,5)


%% EXAMPLE WITH 10 initialStates
rngStream = RandStream('mt19937ar');
for kState = 1:10
    s = rng(kState);
    initState(:,kState) = s.State;
end

rngStream.State = initState(:,1);
rand(rngStream,1,5)
rand(rngStream,1,5)

%reset :: the rand function will output the same numbers as it did at the
%beginning
rngStream.State = initState(:,1);
rand(rngStream,1,5)


%% use an atm object
src = source;
tel = telescope(8,'resolution',128,'samplingTime',0.01);
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0]*1e3,...
    'fractionnalR0',[1],...
    'windSpeed',[10],...
    'windDirection',[0]);

tel = tel + atm;
% display atm
imagesc(tel,src)
% advance phase-screen
for i = 1:10, +tel; imagesc(tel,src), drawnow, end

% set the initialState of your own choosing and restart from the
% initPhaseScreen
while(1)
kState = 10;% user-defined state
atm.initState = initState(:,kState);
reset(tel)
for i = 1:10, +tel; imagesc(tel,src), drawnow, end
end

