clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GECAD Polytechnic of Porto 2020
% Send a mail to Fernando Lezama (ing.flezama@gmail.com) for questions
%This repository contains the experimental set up used for the paper: 
%Bidding in Local Electricity Markets with CascadingWholesale Market Integration
%submitted to IJEPES journal Elsevier 2020

%To replicate the experiments, access to codes folder and run "mainLEM_2020_IJEPES.m" file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;close all;
%tTotalTime=tic; % lets track total computational time
addpath('Functions')
for Select_Algorithm=1:4
    tTotalTime=tic; % lets track total computational time
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load MH parameters (e.g., get MH parameters from DEparameters.m file)
    switch Select_Algorithm
        
        case 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Participants can include their algorithms here
            addpath('HyDE-DF')
            algorithm='VS_EEM'; %'The participants should include their algorithm here'
            DEparameters %Function defined by the participant
            deParameters.I_NP=5; %Notice that some algorithms are limited to one individual
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            deParameters.I_strategy=3; % 1 FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
            % 2: FM_bm=repmat(FVr_bestmemit,I_NP,1);
            % 3: Activates three strategies
            deParameters.I_strategyVersion=1;  % I_strategyVersion==1; %Emulates Vortex algorithm
            % I_strategyVersion==2; %HyDE-DF
            % I_strategyVersion==3; %HyDE
            
        case 2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Participants can include their algorithms here
            addpath('HyDE-DF')
            algorithm='HyDEDF'; %'The participants should include their algorithm here'
            DEparameters %Function defined by the participant
            deParameters.I_NP=5; %Notice that some algorithms are limited to one individual
           
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            deParameters.I_strategy=3; % 1 FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
            % 2 FM_bm=repmat(FVr_bestmemit,I_NP,1);
            %3 : Activates three strategies
            deParameters.I_strategyVersion=2;  % I_strategyVersion==1; %Emulates Vortex algorithm
            % I_strategyVersion==2; %HyDE-DF
            %  I_strategyVersion==3; %HyDE
            
        case 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Participants can include their algorithms here
            addpath('ACO_alg')
            algorithm='ACO_EEM'; %'The participants should include their algorithm here'
            ACOparameters %Function defined by the participant
            acoParameters.maxIter = 500; %in EEM 500
            acoParameters.antNo = 20; %in EEM 20
            acoParameters.rho = 0.3; % inEEM 0.3 Evaporation rate (0.5 standard)
            acoParameters.alpha =0.5; % in EEM 0.5;%.5;  % Phromone exponential parameters
            acoParameters.beta = 0.5;% in EEM 0.5;%0.2;  % Desirability exponetial paramter
            acoParameters.DiscretLevel=30; %in EEM 30
            No_solutions=acoParameters.antNo; %Notice that some algorithms are limited to one individual
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        case 4
             %Participants can include their algorithms here
            addpath('CUMDANCauchyALG')
            algorithm='CUMDANC'; %'The participants should include their algorithm here'
            CUMDANCauchy_parameters
            No_solutions=CUMDANCauchyParameters.I_NP;
        otherwise
            fprintf(1,' No parameters loaded\n');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load Data base
    noRuns=20;
    DB=1; 
    % 1: 9 agents (3 Cons, 3 Pros, 3 Gen) - Classic example
    % 2: 25 consumers (80% PV) and 5 gen - one day
    % 3: 25 consumers (0% PV) and 5 gen - one day
    
    [caseStudyData, DB_name]=callDatabase(DB);
    %% Label of the algorithm and the case study
    Tag.algorithm=algorithm;
    Tag.DB=DB_name;
        
    %% Parameters that can be varied with different inputs
    %EDP comercia tri-horaria
    Ext_supplier=[0.093	0.093	0.093	0.093	0.093	0.093	0.093	0.093	0.16	0.16	0.33	0.33	0.33	0.16	0.16	0.16	0.16	0.16	0.16	0.33	0.33	0.16	0.093	0.093];
    %WS market price MIBEL Average of the Week 5-9 / 08 / 2019
    WS_price=[0.049228	0.047272	0.044936	0.043782	0.043554	0.045226	0.046564	0.048322	0.049872	0.050096	0.049628	0.049796	0.050116	0.050712	0.049648	0.048644	0.048066	0.047712	0.047562	0.048038	0.0491	0.051456	0.050584	0.049478];
    Original=0.25; %0.28; %Grid tariff (Can be modified in function of the case study)
    Feed_in=0.095; %Portugal
    Agg_fee=0.15; %Considered aggregator fee
    
    Gecad_Tariff=[0.1012 0.1012 0.1012 0.1012 0.1012 0.1012 0.1012 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 0.1882 ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Uncomment according to the upper tariff you want to analyze
    %Back_up=repmat(Original,1,caseStudyData.General.numPeriods);
    %Back_up=Ext_supplier;
    Back_up=WS_price+Agg_fee;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Lines used to vary the marginal cost of CHP
    caseStudyData.Type4.MC=caseStudyData.Type4.MC*0;
    caseStudyData.Type4.MC=(caseStudyData.Type4.MC+1)*.1;
    %Line used to increase the load of the case study by a factor
    %caseStudyData.Type1.Load=caseStudyData.Type1.Load.*20; %Increasing the load 20 times
    
       
    %% Set lower/upper bounds of variables
    caseStudyData.General.cf=Feed_in;
    caseStudyData.General.cg=Back_up;
    price_taker=1; %This variable if:
    %0: complete freedom of bidding (FreeM case)
    %1: Consumers with no elasticity of price and quantity (Ptakers case)
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    result_name=['Results_IJEPS_rev1\' algorithm '_' DB_name '_simPV_AGG15WS_PTLM' num2str(price_taker)];
    %result_name=['Results\' algorithm '_' DB_name '_Grid025_PTLM' num2str(price_taker)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lowerB,upperB] = setVariablesBounds(caseStudyData,Back_up,price_taker);
    otherParameters.WCCI_2020_funct='fitnessFun_WCCI2020';
    otherParameters.lowerB=lowerB;
    otherParameters.upperB=upperB;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call the MH for optimizationclear
    ResDB=struc([]);
    for iRuns=1:noRuns %Number of trails
        tOpt=tic;
        rand('state',iRuns)% ensure stochastic indpt trials
        otherParameters.iRuns=iRuns;
        
        switch Select_Algorithm
            case {1, 2}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ResDB(iRuns).Fit_and_p, ...
                    ResDB(iRuns).sol, ...
                    ResDB(iRuns).fitVector]= ...
                    HyDE(deParameters,caseStudyData,otherParameters,lowerB,upperB);
            case 3
                [ResDB(iRuns).Fit_and_p, ...
                    ResDB(iRuns).sol, ...
                    ResDB(iRuns).fitVector,...
                    ResDB(iRuns).learning]= ...
                    ACO_LM(acoParameters,caseStudyData,otherParameters,lowerB,upperB);
             case 4
                [ResDB(iRuns).Fit_and_p, ...
                    ResDB(iRuns).sol, ...
                    ResDB(iRuns).fitVector]= ...
                    CUMDANCauchy_alg(CUMDANCauchyParameters,caseStudyData,otherParameters,lowerB,upperB);
        end
       %[S_val, profits,Penalties, TE]=feval(otherParameters.WCCI_2020_funct,ResDB(iRuns).sol,caseStudyData,otherParameters);
        ResDB(iRuns).tOpt=toc(tOpt); % time of each trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save the results and stats
       
    end
    tTotalTime=toc(tTotalTime); %Total time
    %% End of MH Optimization
    %save(algorithm)
    save(result_name)
    
end





