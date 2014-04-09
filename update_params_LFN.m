function [ params_new ] = update_params_LFN( params, T_count, G_count )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Theta= params.Theta;
Phi = params.Phi;
Beta = params.Beta;
B = params.B;


for p = 1:N
    Theta(p,: ) = T_count(p,:)/sum(T_count(p,:));
end




params_new.Theta =Theta;
params_new.Phi = Phi;
params_new.Beta = Beta;
params_new.B = B;

end

