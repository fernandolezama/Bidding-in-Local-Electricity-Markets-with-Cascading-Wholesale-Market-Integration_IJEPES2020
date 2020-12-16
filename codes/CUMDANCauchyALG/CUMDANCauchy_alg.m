%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:         [FVr_bestmem,S_bestval,I_nfeval] = deopt(fname,S_struct)
% Author:           Rainer Storn, Ken Price, Arnold Neumaier, Jim Van Zandt
% Modified by FLC \GECAD 04/winter/2017
% Modified by TEAM CISECE-UT3/UC-UCLV 2020
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


function [Fit_and_p,solution, fitMaxVector] = ...
    CUMDANCauchy_alg(CUMDANCauchyParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit)

%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = CUMDANCauchyParameters.I_NP;
I_D          = numel(up_habitat_limit); %Number of variables or dimension
I_itermax    = CUMDANCauchyParameters.I_itermax;
selectInd    = CUMDANCauchyParameters.selectInd;

%deParameters.nVariables=I_D;
FVr_minbound = low_habitat_limit;
FVr_maxbound = up_habitat_limit;

%Repair boundary method employed
BRM=CUMDANCauchyParameters.I_bnd_constr; %1: bring the value to bound violated
                               %2: repair in the allowed range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I_strategy   = deParameters.I_strategy; %important variable
fnc= otherParameters.WCCI_2020_funct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----Check input variables---------------------------------------------
if (I_NP < 5)
    I_NP=5;
   fprintf(1,' I_NP increased to minimal value 5\n');
end
%if ((F_CR < 0) || (F_CR > 1))
 %  F_CR=0.5;
  % fprintf(1,'F_CR should be from interval [0,1]; set to default value 0.5\n');
%end
if (I_itermax <= 0)
   I_itermax = 500;
   fprintf(1,'I_itermax should be > 0; set to default value 200\n');
end

%-----Initialize population and some arrays-------------------------------
%FM_pop = zeros(I_NP,I_D); %initialize FM_pop to gain speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-allocation of loop variables
fitMaxVector = nan(1,I_itermax);
% limit iterations by threshold
gen = 1; %iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----FM_pop is a matrix of size I_NPx(I_D+1). It will be initialized------
%----with random values between the min and max values of the-------------
%----parameters-----------------------------------------------------------
% FLC modification - vectorization
minPositionsMatrix=repmat(FVr_minbound,I_NP,1);
maxPositionsMatrix=repmat(FVr_maxbound,I_NP,1);
numParents = CUMDANCauchyParameters.parent;
elitism = CUMDANCauchyParameters.elitism;
%globalBestFitness = inf;
%solution = nan(1, I_NP);

%% generate initial population.
rand('state',otherParameters.iRuns) %Guarantee same initial population
FM_pop=genpop(I_NP,I_D,minPositionsMatrix,maxPositionsMatrix);

[S_val, ~,~]=feval(fnc,FM_pop,caseStudyData,otherParameters);
fitness = sortrows([ (1:I_NP)' S_val(:, 1) ], 2);
%% Selection 
idx = fitness(1:numParents, 1);
FM_parent = FM_pop(idx,:);
globalBestFitness = fitness(1,2);
fitMaxVector(1,gen) = globalBestFitness;
solution = FM_parent(1,:);
%% Store
fprintf('Fitness value: %f\n',fitMaxVector(1,gen) )
fprintf('Generation: %d\n',gen)
gen=gen+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Evaluate the best member after initialization----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------EDA-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------
%FVr_rot  = (0:1:I_NP-1);               % rotating index array (size I_NP)
 
while gen<I_itermax %%&&  fitIterationGap >= threshold
  %------Evaluate the best member after initialization----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Learning
  FM_ui =  learningEDA(FM_parent,I_NP);  
  
%% Boundary Control
  FM_ui = update_eda(FM_ui,minPositionsMatrix,maxPositionsMatrix,BRM); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Evaluation
  [S_val, ~]=feval(fnc,FM_ui,caseStudyData, otherParameters);
  fitness = sortrows([ (1:I_NP)' S_val(:, 1) ], 2);
 
  %% Elitist Selection 
  if elitism > 0 
   idx = fitness(1:elitism, 1);   
   FM_pop(idx,:) = FM_ui(idx,:);
  else 
   FM_pop = FM_ui;
  end
  
  %% parents selection 
  idx = fitness(1:numParents, 1);
  FM_parent = FM_pop(idx,:);
  
  sbestFitness = mean(fitness(1:selectInd,2));
  FM_ui_sel    = FM_ui(fitness(1:selectInd,1),:);
  best         = mean(FM_ui_sel,1);
  
  %if Select_testbed == 1
  % sbest = mean(fitness(1:numParents,2));
  % best  = mean(FM_parent); 
  %else
  % Tesband2
  % sbest = fitness(1,2);
  % best = FM_parent(1,:);
  %end
  
  fitMaxVector(1,gen) =  sbestFitness;
   
   if sbestFitness < globalBestFitness
     globalBestFitness = sbestFitness;
     solution =  best; % best member of current iteration
    end
 
    %% store fitness evolution and obj fun evolution as well
    fprintf('Fitness value: %f %f\n',fitMaxVector(1,gen), globalBestFitness);
    fprintf('Generation: %d\n',gen)

    gen=gen+1;  
    
    %fitMaxVector(1,gen) = sbest;
end %---end while ((I_iter < I_itermax) ...

Fit_and_p=[fitMaxVector(1,gen-1) 0]; %;p2;p3;p4]

solution = update_eda(solution,low_habitat_limit,up_habitat_limit,1);

end

 
% VECTORIZED THE CODE INSTEAD OF USING FOR
function pop=genpop(a,b,lowMatrix,upMatrix)
 pop=unifrnd(lowMatrix,upMatrix,a,b);
end

% VECTORIZED THE CODE INSTEAD OF USING FOR
 function p=update_eda(p,lowMatrix,upMatrix,BRM)
  switch BRM
    case 1 % max and min replace
        [idx] = find(p<lowMatrix);
        p(idx)=lowMatrix(idx);
        [idx] = find(p>upMatrix);
        p(idx)=upMatrix(idx);
    case 2 %Random reinitialization
        [idx] = [find(p<lowMatrix);find(p>upMatrix)];
        replace=unifrnd(lowMatrix(idx),upMatrix(idx),length(idx),1);
        p(idx)=replace;
    case 3 %Bounce Back
        [idx] = find(p<lowMatrix);
      p(idx)=unifrnd(lowMatrix(idx),p(idx),length(idx),1);
        [idx] = find(p>upMatrix);
      p(idx)=unifrnd(p(idx), upMatrix(idx),length(idx),1);
   end
 end

function pop_eda = learningEDA(pop_u,I_NP)
 %Cauchy's distribution
  mu = mean(pop_u);
  sd = std(pop_u);
for i=1:I_NP
 pop_eda(i,:) =  normrnd(mu,sd).*(mu - sd*tan(pi*(rand(1,1))-0.5));
end
end
