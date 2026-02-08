clear
close all
load data2.mat

%% =======================
% USER SETTINGS
%% =======================
nTrainRuns      = 5;   % how many independent training passes
trainingSplits  = [1.0 0.8 0.5 0.3 0.2 0.1];

%% =======================
% DATA PREPROCESSING
%% =======================
T = 5;
fs = 1/T;

Skraplacz = 0.1*data.Skraplacz;
Kompresor1 = data.Kompresor1;
Parownik1 = 0.1*data.Parownik1;
Wlot = 0.1*data.Wlot;
TemperaturaOtoczenia = 0.1*data.TemperaturaOtoczenia;
Rozmraanie = data.Rozmraanie;
WentylatorParownika = data.WentylatorParownika;
ZawrOdszraniania = data.ZawrOdszraniania;

X_all = [ZawrOdszraniania, Kompresor1, WentylatorParownika, Rozmraanie,...
    TemperaturaOtoczenia, Skraplacz, Wlot, Parownik1]';

[X,~] = GetRidOfNans(X_all);
X = X';
n = size(X,2);
t = T*(1:n);

X_toFilter = X(5:8,:);
threshold = [0,5;0,5;-5,0;-5,0];
for i = 1:4
    X_filtered(i,:) = InterpolateData(X_toFilter(i,:),threshold(i,:),t);
end

%% =======================
% SEGMENT DATA INTO CHUNKS
%% =======================
locs = compressorCycles(X(2,:));
locs = locs(2:end);
x = FeatureExtract(X_filtered(3,:),locs);

fridge_air_temp_p2p = x(1,:) - x(2,:);
TF = isoutlier(fridge_air_temp_p2p);

indexContainer = {};
set = 1;
indexContainer{set} = [];

for i = 2:length(x)
    if TF(i)==0
        indexContainer{set} = [indexContainer{set}, locs(i-1):locs(i)];
    else
        if TF(i-1)==0
            set = set + 1;
            indexContainer{set} = [];
        end
    end
end

usableDataIndex = {};
set = 1;
for i = 1:length(indexContainer)
    if length(indexContainer{i}) > 60*2*fs*60
        usableDataIndex{set} = indexContainer{i};
        set = set + 1;
    end
end

nChunks = length(usableDataIndex);




%% =======================
% TRAIN MULTIPLE TIMES
%% =======================
for run = 1:nTrainRuns
    % fprintf('Training run %d / %d\n', run, nTrainRuns);
    %
    % y = cell(1,nChunks);
    % params = cell(1,nChunks);
    %
    % parfor i = 1:nChunks
    %     idx = usableDataIndex{i};
    %     t_sim = T*(0:length(idx)-1);
    %     [y{i}, params{i}] = GATrain( ...
    %         X(2,idx), ...
    %         X_filtered(1,idx), ...
    %         t_sim, ...
    %         X_filtered(3,idx(1)), ...
    %         X_filtered(3,idx));
    % end
    %
    filename = sprintf('training_%d.mat', run);
    % save(filename, 'params', 'y');
    load(filename);
end
results=cell(1,nTrainRuns);


%% =======================
% REPAETE MULTIPLE TIMES DUE TO DATASET SELECTION RANDOMNESS
%% =======================
    AverageResult=cell(1,5);

for recheck=1:5
    ResultMatrix=nan(nTrainRuns,length(trainingSplits));

    for runs = 1:nTrainRuns
        filename = sprintf('training_%d.mat', runs);


        %% =======================
        % EVALUATE TRAINING SPLITS
        %% =======================
        results = struct();

        for s = 1:length(trainingSplits)

            frac = trainingSplits(s);
            nTrain = floor(frac * nChunks);
            nTest  = nChunks - nTrain;

            % Random permutation of chunk indices
            permIdx = randperm(nChunks);

            trainIdx = permIdx(1:nTrain);
            testIdx  = permIdx(nTrain+1:end);

            TrainingData = usableDataIndex(trainIdx);
            TestingData  = usableDataIndex(testIdx);


            fprintf('\nTraining fraction %.0f%% (%d train / %d test)\n', ...
                frac*100, nTrain, nTest);

            %% Load one training run (can average later if desired)
            load(filename, 'params');

            %% === PARAMETER SELECTION (TRAINING ONLY)
            avgRMSE = nan(1,nTrain);

            parfor j = 1:nTrain
                rmse_tmp = nan(1,nTrain);
                for i = 1:nTrain
                    idx = TrainingData{i};
                    t_sim = T*(0:length(idx)-1);
                    rmse_tmp(i) = costFunction( ...
                        -X(2,idx)', X_filtered(1,idx)', ...
                        params{j}(1), params{j}(2), params{j}(3), ...
                        t_sim, X_filtered(3,idx(1)), X_filtered(3,idx));
                end
                avgRMSE(j) = mean(rmse_tmp);
            end

            [~,bestIdx] = min(avgRMSE);
            bestParam = params{bestIdx};

            %% === TESTING (MAPE)

            if length(TestingData) < 1
                TestingData = usableDataIndex;
            end
            pctErr = nan(1,length(TestingData));
            parfor i = 1:length(TestingData)
                idx = TestingData{i};
                t_sim = T*(0:length(idx)-1);

                simData = fridge_fixed_stepCopy( ...
                    X_filtered(1,idx)', ...
                    bestParam(1), bestParam(2), bestParam(3), ...
                    t_sim, X_filtered(3,idx(1))');

                stats = compareOnTime(simData.u1, -X(2,idx));
                pctErr(i) = stats.pctError;
            end

            %% === STORE RESULTS
            results(s).trainingFraction = frac;
            results(s).nTrainingChunks  = nTrain;
            results(s).nTestingChunks   = nTest;
            results(s).nTotalChunks     = nChunks;
            results(s).MAPE             = mean(abs(pctErr));

            fprintf('MAPE = %.3f %%\n', results(s).MAPE);
        end

        %% =======================
        % LEARNING CURVE PLOT
        %% =======================

        trainingFrac = [results.trainingFraction];
        MAPE         = [results.MAPE];
        nTrain       = [results.nTrainingChunks];

        figure
        plot(trainingFrac*100, MAPE, '-o', 'LineWidth', 2)
        grid on

        xlabel('Training data used [%]')
        ylabel('Mean Absolute Percentage Error (MAPE) [%]')
        title('Learning Curve: Training Set Size vs Prediction Error')

        % Optional: annotate number of chunks
        for i = 1:length(nTrain)
            text(trainingFrac(i)*100, MAPE(i), ...
                sprintf('  n=%d', nTrain(i)), ...
                'VerticalAlignment','bottom')
        end

        ResultMatrix(runs,:)=MAPE;
        disp(ResultMatrix);

    end

    AverageResult{recheck}=mean(ResultMatrix,1);
    disp(AverageResult{recheck})
    figure
    plot(trainingSplits*100, AverageResult{recheck}, '-o', 'LineWidth', 2)
    grid on

    xlabel('Training data used [%]')
    ylabel('Mean Absolute Percentage Error (MAPE) [%]')
    title('Learning Curve: Training Set Size vs Prediction Error')

end
% mean(AverageResult)

% Convert cell array to a 10x6 numeric matrix (rows = cells)
M = vertcat(AverageResult{:});

% Compute the mean across the 10 rows to get a 1x6 averaged vector
avgVec = mean(M,1);

% Display result
disp(avgVec);
figure
plot(trainingSplits*100, avgVec, '-o', 'LineWidth', 2)
grid on

xlabel('Training data used [%]')
ylabel('Mean Absolute Percentage Error (MAPE) [%]')
title('Learning Curve: Training Set Size vs Prediction Error')
%%
function simData = fridge_fixed_stepCopy(Ta,G,R,c,t,x0)
% Discrete-time simulation of a refrigerator with hysteresis control
% ––– USER DEFINED SYSTEM ––– ------------------------------------------

% G=2.1686;
% R=48.9120;
% c=114.2959;

A=-1/(R*c);
B=[G/c,1/(R*c)];
C=1;
D=[0,0];

dt     = t(2)-t(1);                % time step [s]
Tend   = t(end);           % simulate 2 hours
% t      = 0:dt:Tend;        % time vector
Ta_fun = @(t) Ta((t/dt)+1);  % outdoor temperature

n   = size(A,1);
x   = zeros(n,length(t));
% x0  = zeros(n,1);         % initial state
x(:,1) = x0;

y  = zeros(1,length(t));
u1 = zeros(1,length(t));
u2 = arrayfun(Ta_fun, t);

% ––– Hysteresis controller state
heater = 0;               % compressor OFF initially
T_low  = -22;               % °C, turn ON below this
T_high = -20;               % °C, turn OFF above this

% ––– MAIN LOOP ––– -----------------------------------------------------
for k = 1:length(t)-1
    y(k) = C*x(:,k) + D*[u1(k); u2(k)];

    % Controller logic (hysteresis)
    if y(k) < T_low
        heater = 0;
    elseif y(k) > T_high
        heater = -1;
    end
    u1(k+1) = heater;

    % Plant update (Euler integration)
    u = [u1(k); u2(k)];
    dx = A*x(:,k) + B*u;
    x(:,k+1) = x(:,k) + dt * dx;
end
y(end) = C*x(:,end) + D*[u1(end); u2(end)];

% ––– OUTPUT ––– --------------------------------------------------------

simData = struct('t',t,'x',x,'y',y,'u1',u1,'u2',u2);

% figure
% subplot(3,1,1)
% plot(t/60,y,'LineWidth',1.5), grid on
% ylabel('Fridge T [°C]')
%
% subplot(3,1,2)
% stairs(t/60,u1,'LineWidth',1.5), grid on
% ylabel('Compressor [0/1]')
%
% subplot(3,1,3)
% plot(t/60,u2,'LineWidth',1.3), grid on
% xlabel('Time [min]')
% ylabel('Outdoor T_a [°C]')
end







