% deParameters.I_NP= 20; % population in DE
% deParameters.F_weight= 0.3; %Mutation factor
% deParameters.F_CR= 0.5; %Recombination constant
% deParameters.I_itermax= 500; % number of max iterations/gen
% deParameters.I_strategy   = 1; %DE strategy
% deParameters.adaptActivated=1; %Change this to one for adaptive DE

deParameters.I_bnd_constr = 3; %Using bound constraints /is possible to change direct in DE
% 1 repair to the lower or upper violated bound 
% 2 rand value in the allowed range
% 3 bounce back

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings used in WCCI paper
deParameters.adaptActivated=0; %Change this to one for adaptive DE
deParameters.I_strategy=3; % 1       FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
                                                    % 2         FM_bm=repmat(FVr_bestmemit,I_NP,1);
                                                    
                                                    
 %3 : Activates three strategies                                                   
deParameters.I_strategyVersion=2;  % I_strategyVersion==1; %Emulates Vortex algorithm
                                                                    % I_strategyVersion==2; %HyDE-DF
                                                                     %  I_strategyVersion==3; %HyDE
deParameters.I_itermax= 2000;
            
deParameters.I_NP=20;
deParameters.F_weight=0.5; %0.5
deParameters.F_CR=0.5; %0.9

