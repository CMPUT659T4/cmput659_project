function loglikelihood = gaussianFit(J, data, lambda)
%% File: gaussianFit.m
%% --------------------------------------------------------------------
%% This function computes the log likelihood of the data given the
%% model, which is defined by its precision matrix J. data is a matrix
%% of N x M, where N is the number of samples (timepoints) and M the 
%% number of variables (Regions of Interest)

n = size(data,1);
emp_cov = (1/n) * data'*data;
loglikelihood = (n/2) * ( log(det(J)) - trace(emp_cov*J) );
