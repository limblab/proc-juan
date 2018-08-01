function [save_dir, data_dir] = get_computer_paths()

[~,hostname]= system('hostname');

switch strtrim(hostname)
    case 'FSMC17SW0D9GTF1' % my 2016 Macbook Pro
        save_dir = '/Users/mattperich/Dropbox/Research/Papers/Juan and Matt - Stability latent activity/Results/';
        data_dir = '/Users/mattperich/Dropbox/Research/Data/TrialDataFiles/';
        
    otherwise
        save_dir = '/Users/juangallego/Dropbox/Juan and Matt - Stability latent activity/Results/';
        data_dir = '/Users/juangallego/Documents/NeuroPlast/Data/';
end

