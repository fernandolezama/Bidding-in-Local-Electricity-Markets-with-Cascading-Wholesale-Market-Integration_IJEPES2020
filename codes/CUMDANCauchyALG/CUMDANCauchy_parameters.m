
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GECAD - GECCO and CEC 2020 Competition: Evolutionary Computation in Uncertain Environments: A Smart Grid Application 
%
% TEAM CISECE-UT3/UC-UCLV
% CUMDANCauchy++: a Cellular EDA Designed to Solve the Energy Resource Management Problem Under Uncertainty
%
% PhD         Ansel Y. Rodriguez-Gonzalez,         ansel@cicese.mx             [1]
% MS  Student Bryan Rodrigo Quiroz Palominos       quirozb@cicese.edu.mx       [1]
% MS  Student Carlos A. Oliva Moreno               carlosalberto@cicese.edu.mx [1]
% Phd Student Yoan Martinez Lopez,                 yoan.martinez@reduc.edu.cu  [2,3]
% PhD         Julio Madera,                        julio.madera@reduc.edu.cu   [2]
%
% [1] CISECE-UT3 (Unidad de Transferencia Tecnológica Tepic del Centro de Investigación Científica y de Educación Superior de Ensenada, Mexico)
% [2] UC (Universidad de Camaguey, Cuba)
% [3] UCLV (Universidad de Central de las Villas, Cuba)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Testband 1;
% CUMDANCauchyParameters.I_NP= 150; % population in DE 10
% CUMDANCauchyParameters.parent = 16;
% CUMDANCauchyParameters.elitism = 0;
% CUMDANCauchyParameters.I_itermax =  floor((50000)/(CUMDANCauchyParameters.I_NP*10))-1;
% CUMDANCauchyParameters.selectInd = 16;

%Testband 2;
CUMDANCauchyParameters.I_NP= 150; % population in DE 10
CUMDANCauchyParameters.parent = 16;
CUMDANCauchyParameters.elitism = 1;
CUMDANCauchyParameters.I_itermax =  floor((10000)/(CUMDANCauchyParameters.I_NP))-1;
CUMDANCauchyParameters.selectInd = 1; 


CUMDANCauchyParameters.I_bnd_constr = 1; %Using bound constraints 
% 1 repair to the lower or upper violated bound 
% 2 rand value in the allowed range
% 3 bounce back


