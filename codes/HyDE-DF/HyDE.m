%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:          Fernando Lezama
% Modified by FLC \GECAD 2020

%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB source code of HyDE-DF which is an improved version of DE with self-adaptive parameters.
%% About HyDE-DF, please see following papers:
%%
%% * Lezama, F., Faia, R., Soares, J., and Vale, Z.: HyDE-DF: A novel self-adaptive version of differential1evolution for numerical optimization
%% * Lezama, F., Soares, J., Faia, R., Pinto, T., and Vale, Z.: A New Hybrid-Adaptive Differential Evolution for a Smart Grid Application Under Uncertainty. In IEEE Congress on Evolutionary Computation (WCCI IEEE-CEC), Rio de Janeiro, Brazil, July 2018.
%% * Lezama, F., Soares, J., Faia, R., and Vale, Z.:Hybrid-Adaptive Differential Evolution with Decay Function(HyDE-DF) Applied to the 100-Digit Challenge Competition onSingle Objective Numerical Optimization

%% For this package, we downloaded DE's source code from Rainer Storn
% (http://www.icsi.berkeley.edu/~storn/code.html)and modified it.
% 
%%
%%%%%%%%%%%%%%%%%%% 

function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
HyDE(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit)


 % HyDE(deParameters,otherParameters,low_habitat_limit,up_habitat_limit,initialSolution)

%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP;
F_weight     = deParameters.F_weight;
F_CR         = deParameters.F_CR;
I_D          = numel(up_habitat_limit); %Number of variables or dimension
deParameters.nVariables=I_D;
FVr_minbound = low_habitat_limit;
FVr_maxbound = up_habitat_limit;
I_itermax    = deParameters.I_itermax;

%Repair boundary method employed
BRM=deParameters.I_bnd_constr; %1: bring the value to bound violated
                               %2: repair in the allowed range
                               %3: Bounce-back

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_strategy   = deParameters.I_strategy; %important variable
fnc= otherParameters.WCCI_2020_funct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----Check input variables---------------------------------------------
if (I_NP < 5)
   I_NP=5;
   fprintf(1,' I_NP increased to minimal value 5\n');
end
if ((F_CR < 0) || (F_CR > 1))
   F_CR=0.5;
   fprintf(1,'F_CR should be from interval [0,1]; set to default value 0.5\n');
end
if (I_itermax <= 0)
   I_itermax = 200;
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
deParameters.minPositionsMatrix=minPositionsMatrix;
deParameters.maxPositionsMatrix=maxPositionsMatrix;
% generate initial population.
rand('state',otherParameters.iRuns) %Guarantee same initial population
FM_pop=genpop(I_NP,I_D,minPositionsMatrix,maxPositionsMatrix);

%% Correct the solution to avoid selling to a cheaper price than the marginal cost
T=caseStudyData.General.numPeriods; %Periods
numType1=caseStudyData.General.numN1; %Consumers only
numType2=caseStudyData.General.numN2; %Prosumers

ind_price=numel(up_habitat_limit)/2;
Q_lead=-1*FM_pop(:,T*(numType1+numType2)+1:ind_price);
if isempty(Q_lead)~=1
    b_chp=caseStudyData.Type4.MC; % cost factor for the CHP generation
    c_CHP=@(b_chp,g_chp) ((b_chp.*sqrt(g_chp))./g_chp+1e-10);
    Marginal_Cost=c_CHP(b_chp(1,1),Q_lead);

    minPossibleBid=minPositionsMatrix(:,ind_price+1+T*(numType1+numType2):end);
    New_bid=max(Marginal_Cost,minPossibleBid);
    Corrected_bid= unifrnd(New_bid,maxPositionsMatrix(:,ind_price+1+T*(numType1+numType2):end))   ;

    ind_to_zero=isnan(Corrected_bid);
    Q_lead(ind_to_zero)=0;
    Corrected_bid(ind_to_zero)=caseStudyData.General.cf;

    FM_pop(:,T*(numType1+numType2)+1:ind_price)=-1*Q_lead;
    FM_pop(:,ind_price+1+T*(numType1+numType2):end)=Corrected_bid;
end


if nargin>5
    noInitialSolutions = size(initialSolution,1);
    FM_pop(1:noInitialSolutions,:)=initialSolution;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Evaluate the best member after initialization----------------------
% Modified by FLC
[S_val, ~]=feval(fnc,FM_pop,caseStudyData,otherParameters);

%S_val=feval(fnc,FM_pop);

[~,I_best_index] = min(S_val); % This mean that the best individual correspond to the best worst performance
FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration

fitMaxVector(:,1)=S_val(I_best_index); %We save the mean value and mean penalty value
% The user can decide to save the mean, best, or any other value here

%------DE-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------
FVr_rot  = (0:1:I_NP-1);               % rotating index array (size I_NP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYDE
if deParameters.I_strategy==3
    F_weight_old=repmat(F_weight,I_NP,3);
    F_weight= F_weight_old;
    F_CR_old=repmat(F_CR,I_NP,1);
    F_CR=F_CR_old;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_strategyVersion=deParameters.I_strategyVersion;


while gen<=I_itermax  %%&&  fitIterationGap >= threshold
    %a = itr / MaxItr; % a value for gammaincinv function
    other.a=(I_itermax-gen)/I_itermax;
    other.lowerlimit=minPositionsMatrix; %lower limit of the problem (in matrix form)
    other.upperlimit = maxPositionsMatrix; %upper limit of the problem (in matrix form)
    
    
     if deParameters.I_strategy==3 %HyDE and HyDE-DF (jDE adaptive mechanism)
              
                value_R=rand(I_NP,3);
                ind1=value_R<0.1;
                ind2=rand(I_NP,1)<0.1;
                F_weight(ind1)=0.1+rand(sum(sum(ind1)),1)*0.9;
                F_weight(~ind1)=F_weight_old(~ind1);
                F_CR(ind2)=rand(sum(ind2),1);
                F_CR(~ind2)=F_CR_old(~ind2);
            
     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [FM_ui,FM_base,~]=generate_trial(I_strategy,F_weight, F_CR, FM_pop, FVr_bestmemit,I_NP, I_D, FVr_rot,I_strategyVersion,other);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

    %% Boundary Control
    FM_ui=update(FM_ui,minPositionsMatrix,maxPositionsMatrix,BRM,FM_base);
       %% Correct the solution to avoid selling to a cheaper price than the marginal cost
    Q_lead=-1*FM_ui(:,T*(numType1+numType2)+1:ind_price);
    if isempty(Q_lead)~=1
        Marginal_Cost=c_CHP(b_chp(1,1),Q_lead);

        minPossibleBid=minPositionsMatrix(:,ind_price+1+T*(numType1+numType2):end);
        New_bid=max(Marginal_Cost,minPossibleBid);
        Corrected_bid= unifrnd(New_bid,maxPositionsMatrix(:,ind_price+1+T*(numType1+numType2):end))   ;    ind_to_zero=isnan(Corrected_bid);

        Q_lead(ind_to_zero)=0;
        Corrected_bid(ind_to_zero)=caseStudyData.General.cf;
        FM_ui(:,T*(numType1+numType2)+1:ind_price)=-1*Q_lead;
        FM_ui(:,ind_price+1+T*(numType1+numType2):end)=Corrected_bid;
    end
%     
    

    %%Evaluation of new Pop
    %S_val_temp=feval(fnc,FM_ui);
    [S_val_temp, ~]=feval(fnc,FM_ui,caseStudyData,otherParameters);
    
    
    %% Elitist Selection
    ind=find(S_val_temp<S_val);
    S_val(ind)=S_val_temp(ind);
    FM_pop(ind,:)=FM_ui(ind,:);
  
  
    %% update best results
    [S_bestval,I_best_index] = min(S_val);
    FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
    % store fitness evolution and obj fun evolution as well
    fitMaxVector(1,gen) = S_bestval;
    %S_bestval
    
     if deParameters.I_strategy==3 %jDE mechanism
        F_weight_old(ind,:)=F_weight(ind,:);
        F_CR_old(ind)=F_CR(ind);
     end
    
     
      fprintf('Fitness value: %f\n',fitMaxVector(1,gen) )
    fprintf('Generation: %d\n',gen)
   
    gen=gen+1;
    %S_bestval
    
    

end %---end while ((I_iter < I_itermax) ...
%Fit_and_p=fitMaxVector(1,gen-1); 
%p1=0;
Fit_and_p=[fitMaxVector(1,gen-1) 0]; %;p2;p3;p4]

 
% VECTORIZED THE CODE INSTEAD OF USING FOR
function pop=genpop(a,b,lowMatrix,upMatrix)
pop=unifrnd(lowMatrix,upMatrix,a,b);

% VECTORIZED THE CODE INSTEAD OF USING FOR
function p=update(p,lowMatrix,upMatrix,BRM,FM_base)
switch BRM
    case 1 %Our method
        %[popsize,dim]=size(p);
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
      p(idx)=unifrnd(lowMatrix(idx),FM_base(idx),length(idx),1);
        [idx] = find(p>upMatrix);
      p(idx)=unifrnd(FM_base(idx), upMatrix(idx),length(idx),1);
end

function [FM_ui,FM_base,msg]=generate_trial(method,F_weight, F_CR, FM_pop, FVr_bestmemit,I_NP,I_D,FVr_rot,I_strategyVersion,other)

  FM_popold = FM_pop;                  % save the old population
  FVr_ind = randperm(4);               % index pointer array
  FVr_a1  = randperm(I_NP);                   % shuffle locations of vectors
  FVr_rt  = rem(FVr_rot+FVr_ind(1),I_NP);     % rotate indices by ind(1) positions
  FVr_a2  = FVr_a1(FVr_rt+1);                 % rotate vector locations
  FVr_rt  = rem(FVr_rot+FVr_ind(2),I_NP);
  FVr_a3  = FVr_a2(FVr_rt+1);                
  FM_pm1 = FM_popold(FVr_a1,:);             % shuffled population 1
  FM_pm2 = FM_popold(FVr_a2,:);             % shuffled population 2
  FM_pm3 = FM_popold(FVr_a3,:);             % shuffled population 3
  %FM_mui = rand(I_NP,I_D) < F_CR;  % all random numbers < F_CR are 1, 0 otherwise
  %FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
  
    if length(F_CR)==1  %Meaning the same F_CR for all individuals
        FM_mui = rand(I_NP,I_D) < F_CR;  % all random numbers < F_CR are 1, 0 otherwise
        FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
    else %Meaning a different F_CR for each individual
        FM_mui = rand(I_NP,I_D) < repmat(F_CR,1,I_D);  % all random numbers < F_CR are 1, 0 otherwise
        FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
    end


    switch method
        case 1
            FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;     % crossover
            FM_base = FM_pm3;
            msg=' DE/rand/bin';
        case 2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            %VEC by FLC
            FM_bm=repmat(FVr_bestmemit,I_NP,1);
            FM_ui = FM_popold + F_weight*(FM_bm-FM_popold) + F_weight*(FM_pm1 - FM_pm2);
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
            FM_base = FM_bm;
            msg=' DE/current-to-best/1';
        
        
        case 3 %jDEPerturbated_v3 v4... v7
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            FM_bm=repmat(FVr_bestmemit,I_NP,1);
            
            if  I_strategyVersion==1 %Emulates Vortex algorithm
                        a=other.a;
                        %a = itr / MaxItr; % a value for gammaincinv function
                        %a_par(count)=a;
                        ginv = (1/0.1)*gammaincinv(0.1,a); % compute the new ginv value
                        %ginv_par(count)=ginv;
                        r = ginv * ((other.upperlimit - other.lowerlimit) / 2); %decrease the radius
                        C = r.*randn(I_NP,I_D);
                        FM_ui = bsxfun(@plus, C, FM_bm(1,:));
            end
            
            if  I_strategyVersion==2 %HyDE-DF
                a=other.a; %Linear decrease
                ginv=exp((1-(1/a^2))); %Exponential decreasing funtion
                %Not so fast but convergence
                FM_ui = FM_popold + repmat(F_weight(:,3),1,I_D).*(FM_pm1 - FM_pm2)  + ginv*(repmat(F_weight(:,1),1,I_D).*(FM_bm.*(repmat(F_weight(:,2),1,I_D)+randn(I_NP,I_D))-FM_popold));   % differential variation
            end
                
            if  I_strategyVersion==3 %HyDE
                FM_ui = FM_popold + repmat(F_weight(:,1),1,I_D).*(FM_bm.*(repmat(F_weight(:,2),1,I_D)+randn(I_NP,I_D))-FM_popold) + repmat(F_weight(:,3),1,I_D).*(FM_pm1 - FM_pm2);
            end
              
        
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
            FM_base = FM_bm;
            msg=' HyDE/current-to-best/1';   
            
    end
return

