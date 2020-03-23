function [c_E_hat,Em,sse] = EIOT_PRED(S_E,NUM_NC,dm,c_A_bounds,beq)
%Inputs:
% S_E : Extended IOT spectral vector                        [lambda x n ]
% NUM_NC: number of non-chemical interference component,it has to be a scalar value,
% the last x number of column of S_E is the non-chemical interference. 
% You can set this to zero if you don't have any vectors for non-chemical interference.
% dm  : Spectral measurement                                [ lambda x 1 ]
% c_A_bounds: Min and Max of allowed concentration for those non-chemical interference vectors. 
% initially set to [-inf inf]. You can set this to empty if no non-chemical interference is needed.
%Outputs:
% c_E_hat: where the first n elements represent the mass fraction of
% chemical species in mixture, element beyond n is the strength for non-chemical
% interference
% Em: original residual spectrum
% sse : Squared Spectral Error, sum of squares of residuals after EIOT
% deflation.
% 
% Created by Sal Garcia
% Modified by Zhenqi (Pete) Shi
% Eli Lilly and Company
% 2018.8.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


H   =  S_E'*S_E;
f   = -S_E'* dm;
A   = [];
b   = [];
if length(c_A_bounds)~=0
    Aeq = [ones(1,size(S_E,2)-NUM_NC) zeros(1,NUM_NC)];
    % beq = 1;
    lb = [zeros(1,size(S_E,2)-NUM_NC),ones(1,NUM_NC)*c_A_bounds(1)];
    ub = [ones(1,size(S_E,2)-NUM_NC) ,ones(1,NUM_NC)*c_A_bounds(2)];
else
    Aeq = [ones(1,size(S_E,2))];
    % beq = 1;
    lb = [zeros(1,size(S_E,2))];
    ub = [ones(1,size(S_E,2)) ];
end


c_E_hat = quadprog(H,f,A,b,Aeq,beq,lb,ub);
dm_hat  = S_E*c_E_hat;
Em      = dm-dm_hat;
sse     = sum((dm - dm_hat).^2);
