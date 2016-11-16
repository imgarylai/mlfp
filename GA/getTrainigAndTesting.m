function [train,test] = getTrainigAndTesting(d) 
% Split data into training set and testing set
% Select 'perfect' samples as testing set where the number of NaN is less than 4
% Output will be a struct that contains everything

d_copy=d;

testing_keep_indices = (sum(isnan(d.Features), 1) < 4);
training_keep_indices = ~testing_keep_indices;

d.Features = d.Features(:,testing_keep_indices); 
d.Survival = d.Survival(testing_keep_indices);
d.Censored = d.Censored(testing_keep_indices);
d.AvailableProtein = d.AvailableProtein(testing_keep_indices);
d.AvailableClinical = d.AvailableClinical(testing_keep_indices);
d.AvailableCNV = d.AvailableCNV(testing_keep_indices);
d.AvailableMutation = d.AvailableMutation(testing_keep_indices);
d.AvailablemRNA = d.AvailablemRNA(testing_keep_indices);
test = d;

d_copy.Features = d_copy.Features(:,training_keep_indices); 
d_copy.Survival = d_copy.Survival(training_keep_indices);
d_copy.Censored = d_copy.Censored(training_keep_indices);
d_copy.AvailableProtein = d_copy.AvailableProtein(training_keep_indices);
d_copy.AvailableClinical = d_copy.AvailableClinical(training_keep_indices);
d_copy.AvailableCNV = d_copy.AvailableCNV(training_keep_indices);
d_copy.AvailableMutation = d_copy.AvailableMutation(training_keep_indices);
d_copy.AvailablemRNA = d_copy.AvailablemRNA(training_keep_indices);

train = d_copy;