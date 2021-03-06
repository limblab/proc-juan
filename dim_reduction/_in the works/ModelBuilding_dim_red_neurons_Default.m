function options = ModelBuilding_dim_red_neurons_Default()
%       options             : structure with fields
%           fillen              : filter length in seconds (tipically 0.5)
%           UseAllInputs        : 1 to use all inputs, 0 to specify a neuronID file, or a NeuronIDs array
%           PolynomialOrder     : order of the Weiner non-linearity (0=no Polynomial)
%           PredEMG, PredForce, PredCursPos, PredVeloc, numPCs :
%                                 flags to include EMG, Force, Cursor Position
%                                 and Velocity in the prediction model
%                                 (0=no,1=yes), if numPCs is present, will
%                                 use numPCs components as inputs instead of
%                                 spikeratedata
%           Use_Thresh,Use_EMGs,Use_Ridge,Use_SD:
%                                 options to fit only data above a certain
%                                 threshold, use EMGs as inputs instead of
%                                 spikes, or use a ridge regression to fit model
%                                 
%           EMGcascade          : in mfxval, will call BuildModel twice, to
%                                 provide a neuron-to-emg decoder followed 
%                                 by an emg-to-cursorposition decoder
%           plotflag            : plot predictions after xval

options = struct(...
    'nbr_dims', 2,...
    'PredEMGs', 0,...
    'PredForce',0,...
    'PredCursPos',0,...
    'PredVeloc',0,...
    'fillen',0.5,...
    'UseAllInputs',1,...
    'PolynomialOrder',0,...
    'numPCs',0,...
    'Use_Thresh',0,...
    'Use_EMGs',0,...
    'Use_Ridge',0,...
    'Use_SD',0,...
    'EMGcascade',0,...
    'plotflag',0,...
    'foldlength',60);
