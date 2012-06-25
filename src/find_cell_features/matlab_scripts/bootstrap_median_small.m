function range = bootstrap_mean_small(data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('data',@isnumeric);

i_p.addParamValue('bonfer_correction',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('samp_num',100000,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('samp_size',50,@(x)isnumeric(x) && x > 0);

i_p.parse(data,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samp_num = i_p.Results.samp_num;
samp_size = i_p.Results.samp_size;

% sample_indexes = NaN(samp_num,samp_size);
% for i=1:samp_num
%     sample_indexes(i,:) = randi(length(data),[samp_size,1]);
% end

% tic;
% if (mod(samp_size,1) == 0)
%     median_index = [samp_size/2,samp_size/2+1];
% else
%     median_index = samp_size/2+0.5;
% end
% 
% data_sort = sort(data);
% samples = NaN(samp_num,1);
% for i=1:samp_num
%     temp = sort(randi(length(data),[samp_size,1]));
%     temp = temp(median_index);
%     samples(i) = median(data_sort(temp));
% end
% toc;

% tic;
samples = NaN(samp_num,1);
for i=1:samp_num
    samples(i) = median(data(randi(length(data),[samp_size,1])));
end
% toc;

% round(length(samples)*0.05/i_p.Results.bonfer_correction)

samples = sort(samples);

range = [samples(round(length(samples)*0.05/i_p.Results.bonfer_correction)), ... 
    samples(round(length(samples)*(1-0.05/i_p.Results.bonfer_correction)))];