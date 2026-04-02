% Code modified from that used in Vidaurre et al. (2017) PNAS

%% setup the MATLAB path
addpath(genpath(['path/to/HMM-MAR']))

% Load the filenames and paths
H2 = importdata('path/to/filenames.txt');

% Pre-define constants
K = 12;  % Number of states
repetitions = 1;  % Number of repetitions
DirData = 'path/to/in';
DirOut = 'path/to/out';
TR = 0.735;
use_stochastic = 1;

% Number of subjects
N = ?;

% Preallocate memory for cell arrays
f = cell(N,1); 
T = cell(N,1);

% Enable parallel pool if not already active
if isempty(gcp('nocreate'))
    parpool;
end

parfor j = 1:N
    % Open the file
    fid = fopen([H2{j} '/dr_stage1.txt'], 'r');
    
    % Check if file opened successfully
    if fid == -1
        error('File could not be opened.');
    end
    
    % Read the file content with fscanf (formatted reading)
    % Adjust the format specifier based on your file structure
    data = fscanf(fid, '%f', [100, inf]);  % Specify the number of columns in the file
    
    % Close the file
    fclose(fid);
    
    % Process the data (first 490 columns)
    f{j} = data(:, 1:490)';
    T{j} = 490;  % Set the time points length
    
    % Select the columns you are interested in
    f{j} = f{j}(:,[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 45 46 48 49 50 52 53 57 58 60 63 64 93]);
end

options = struct();
options.K = K; % number of states 
options.order = 0; % no autoregressive components
options.zeromean = 0; % model the mean
options.covtype = 'full'; % full covariance matrix
options.Fs = 1/TR; 
options.verbose = 1;
options.standardise = 1;
options.inittype = 'HMM-MAR';
options.cyc = 300;
options.initcyc = 10;
options.initrep = 3;
options.dropstates =1;
options.useParallel = 1;

% stochastic options
if use_stochastic
    options.BIGNbatch = round(N/100);
    options.BIGtol = 1e-5;
    options.BIGcyc = 600;
    options.BIGundertol_tostop = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

% We run the HMM multiple times
for r = 1:repetitions
 [hmm, Gamma, ~, vpath] = hmmmar(f,T,options);
 save([DirOut 'HMMrun_rep' num2str(r) '.mat'],'Gamma','vpath','hmm')
 disp(['RUN ' num2str(r)])
end

save('biobank_hmm.mat')
