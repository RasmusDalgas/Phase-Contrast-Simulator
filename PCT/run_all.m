close all; clear; clc;

% run all simulations and save results in folder (can take hours)

for cc = 5:6
    
  eval(sprintf('compare_TSD_ACD0%d',cc))
    
end

