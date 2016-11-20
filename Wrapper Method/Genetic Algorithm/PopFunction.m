function [pop] = PopFunction(GenomeLength,~,options)
% Generate selected feature index as a logical vector 
% For exmaple, [0 1 0 0 1 ..... 1 1 0 0]
rand('seed',1);
RD = rand;  
pop = (rand(options.PopulationSize, GenomeLength)> RD); % Initial Population
end
