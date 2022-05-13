function enableFieldTrip(type,path)
% Adds or removes Fieldtrip toolbox
%
% Scripts is based on online recommendations from
% 'http://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path'

% path = 'C:\Users\wurmf\Documents\Github\fieldtrip\';

if type == 1    
    addpath(path)   
    ft_defaults   
elseif type == 0
    rmpath(genpath(path))
else
    warning('no type specified!!! 1 - add fieldtrip, 0 - remove fieldtrip')
end

