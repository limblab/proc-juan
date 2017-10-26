%
% Test any model on any dataset
%
% This is a quick hack with many things to do: 
%   1) make it compatible with Matt's evalModel and getModel
%   2) Allow for different evaluation metrics, xval...

function [r2, trial_data] = testModel(trial_data,model_info)


% Rename params
in_signals = model_info.in_signals;
out_signals = model_info.out_signals;
b = model_info.b;
model_name = model_info.model_name;
model_type = model_info.model_type;


% Add predictions to trial_data
for trial = 1:length(trial_data)
    
    x  = get_vars(trial_data(trial),in_signals);
    yfit = zeros(size(x,1),size(b,2));
    
    for iVar = 1:size(b,2)
        % the model B is built so as the first term is the offset
        yfit(:,iVar) = [ones(size(x,1),1), x]*b(:,iVar);
    end
    
    trial_data(trial).([model_type '_' model_name]) = yfit;
end


% Compute R2
y_test = get_vars(trial_data,out_signals);
y_fits = get_vars(trial_data,{[model_type '_' model_name],out_signals{2}});
r2 = compute_r2(y_test,y_fits);
% r2 = calc_r(y_test,y_fits);