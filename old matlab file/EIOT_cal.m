function OUTPUT=EIOT_CAL(X_CAL,Y_CAL,K_AUG,X_VAL,Y_VAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original CLS and PACLS code
% Inputs
%   X_CAL and Y_CAL are training data for spectra and concentration
%   K_AUG: augmentation vector if any
%   X_VAL and Y_VAL are test data for spectra and concentration. If only
%   X_CAL and Y_CAL are available, please feel free to use X_CAL and Y_CAL
%   for X_VAL and Y_VAL.
% Outputs
%   K: estimated pure component spectra for chemical species or it is a
%   concatenation of K and K_AUG (see line 25)
%   All relevant outputs, see bellows.
%
% Zhenqi (Pete) Shi
% Duquesne University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%K=inv(Y_CAL'*Y_CAL)*Y_CAL'*X_CAL;
K=pinv(Y_CAL)*X_CAL;
size(K);
Y_CAL_PRED=X_CAL*pinv(K);

RES_CAL=(Y_CAL*K-X_CAL);

K_FIN=[K;K_AUG];

Y_VAL_PRED=X_VAL*pinv(K_FIN);

RES_VAL=(Y_VAL_PRED*K_FIN-X_VAL);

for i=1:size(Y_CAL,2)
    RMSEC(i)=sqrt(sum((Y_CAL(:,i)-Y_CAL_PRED(:,i)).^2)/size(Y_CAL,1));
    RMSEP(i)=sqrt(sum((Y_VAL(:,i)-Y_VAL_PRED(:,i)).^2)/size(Y_VAL,1));
end

OUTPUT.K=K;
OUTPUT.Y_CAL_PRED=Y_CAL_PRED;
OUTPUT.Y_VAL_PRED=Y_VAL_PRED;
OUTPUT.RMSEC=RMSEC;
OUTPUT.RMSEP=RMSEP;
OUTPUT.RES_CAL=RES_CAL;
OUTPUT.RES_VAL=RES_VAL;

