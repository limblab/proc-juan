function [save_dir, data_dir] = get_computer_paths()

[~,hostname]= system('hostname');

switch strtrim(hostname)
    case 'FSMC17SW0D9GTF1' % my 2016 Macbook Pro
        save_dir = '/Users/mattperich/Dropbox/Research/Papers/Juan and Matt - Stability latent activity/Results/';
        data_dir = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/';
        
    otherwise % assume it's Juan's until he fills in his computer
        save_dir = '/Users/juangallego/Dropbox/M1 FES Grant documents/Juan 2018 - Stable latent activity/Results/';
        data_dir = '/Users/juangallego/Documents/NeuroPlast/Data/TrialDataFiles/';
end

