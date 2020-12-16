clc;clear;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GECAD Polytechnic of Porto 2020
% Send a mail to Fernando Lezama (ing.flezama@gmail.com) for questions
%This repository contains the experimental set up used for the paper: 
%Bidding in Local Electricity Markets with CascadingWholesale Market Integration
%submitted to IJEPES journal Elsevier 2020

%To replicate the experiments, access to codes folder and run "mainLEM_2020_IJEPES.m" file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tTotalTime=tic; % lets track total computational time
%addpath('CallDataBases')
No_Algorithm=4;
%DB_name='9agents';
DB_name='25C_0P_5CHP_1day'
%DB_name='5C_20P_5CHP_1day'
PT=1;
DB=3; %1 (500) and 2 (1); %Select the database you want to analyze

for Select_Algorithm=1:No_Algorithm
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load MH parameters (e.g., get MH parameters from DEparameters.m file)
    switch Select_Algorithm
        case 1
            alg='ACO_EEM'
        case 2
            alg='HyDEDF'   
        case 3
            alg='VS_EEM'
        case 4
            alg='CUMDANC'
            
        otherwise
            fprintf(1,' No parameters loaded\n');
    end
     algorithm=['Results_IJEPS_rev1\' alg '_' DB_name '_AGG15WS_PTLM' num2str(PT)]; %'The participants should include their algorithm here'
    
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load Data base
    noRuns=20;
    %DB=1; %1 (500) and 2 (1); %Select the database you want to analyze
    % 1: playing example
    % 2: 10 agents real data
    % 3: 100 agents real data
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
    price_taker=PT; %This variable if 0, normal consumers, if 1: consumers with no elasticity of demand and price in the LM
    %0: complete freedom of bidding
    %1: Consumers with no elasticity of price and quantity (max buy price)
    %2: prosumers with no elasticity of price and quantity (min sell price)
    [lowerB,upperB] = setVariablesBounds(caseStudyData,Back_up,price_taker);
    %% With this the LM is automatically dissable
    lowerB=lowerB*0;
    upperB=upperB*0;
    otherParameters.WCCI_2020_funct='fitnessFun_WCCI2020';
    otherParameters.lowerB=lowerB;
    otherParameters.upperB=upperB;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call the MH for optimizationclear
    ResDB=struc([]);
    %% Making nice pictures for the IJEPS2020
    load(algorithm, 'ResDB')
    
%% Convergency I cannot use this because different number of iterations
%     for iRuns=1:noRuns %Number of trails
%         val=ResDB(iRuns).fitVector;
%         Matrix_convergency(iRuns,:,Select_Algorithm)=val;
%     end
    %find the best result
    best_val=Inf;
    for iRuns=1:noRuns %Number of trails
        val=ResDB(iRuns).Fit_and_p;
        if val<best_val
            best_ind=iRuns;
            best_val=val;
        end
    end
    sol_test=ResDB(best_ind).sol; %We can take the best solution found
    %sol_test=X
    [S_val_best, profits_best,Penalties_best, TE_best, LM_CP] = feval('fitnessFun_WCCI2020',sol_test,caseStudyData,otherParameters);
    S_val_best
    T_profits=sum(profits_best)
    MH_individual_prof(:,Select_Algorithm)=sum(TE_best.profits,2);
    
    
    for iRuns=1:noRuns %Number of trails
        sol_test=ResDB(iRuns).sol; %We can take the best solution found
        [S_val, profits,Penalties, TE] = feval('fitnessFun_WCCI2020',sol_test,caseStudyData,otherParameters);
        ResDB(iRuns).Fit_and_p
        S_val
        profits_a=sum(TE.profits,2);
        M_profitT(iRuns)=sum(profits); %Total profits of the system
        M_fitness(iRuns)=S_val; %Fitness of the system
        M_profitCons(iRuns)=sum(profits_a(1:25,1));
        M_profitProsu(iRuns)=sum(profits_a(26:26,1));
        M_profitProd(iRuns)=sum(profits_a(27:31,1));
        M_Time(iRuns)=ResDB(iRuns).tOpt;
    end
    Table(:,:,Select_Algorithm)=[ M_profitT', M_profitCons' ,M_profitProsu',M_profitProd',M_fitness',M_Time'/60 ];
    
    %Recheck best solution due to error in CUDAMcauchy store system
    [best_val,best_ind]=min(M_fitness);
    sol_test=ResDB(best_ind).sol; %We can take the best solution found
    %sol_test=X
    [S_val_best, profits_best,Penalties_best, TE_best, LM_CP] = feval('fitnessFun_WCCI2020',sol_test,caseStudyData,otherParameters);
    S_val_best
    T_profits=sum(profits_best)
    MH_individual_prof(:,Select_Algorithm)=sum(TE_best.profits,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %This information is saved for the best value only
    Table_CPrices(Select_Algorithm,:)=LM_CP;
    
    for i=1:length(profits) % Calculate the energy transacted in the LM and in the WS/retail market
        energy_LM(:,i)=TE_best.TransE(:,1,i); %[LM; GridT; demand/generation]
        energy_Grid(:,i)=TE_best.TransE(:,2,i); %[LM; GridT; demand/generation]
    end
    CON_mask_LM=energy_LM<0; %Buy
    CON_mask_Grid=energy_Grid<0; %Buy
    
    GEN_mask_LM=energy_LM>0; %Sell
    GEN_mask_Grid=energy_Grid>0; %Sell
        
    LM_CON_B=sum(energy_LM.*CON_mask_LM,2);
    WS_CON_B=sum(energy_Grid.*CON_mask_Grid,2);
    LM_GEN_S=sum(energy_LM.*GEN_mask_LM,2);
    WS_GEN_S=sum(energy_Grid.*GEN_mask_Grid,2);
    
    Matrix_energyBuy(:,:,Select_Algorithm)=[LM_CON_B';WS_CON_B'];
    Matrix_energySell(:,:,Select_Algorithm)=[LM_GEN_S';WS_GEN_S'];
    
end

for Select_Algorithm=1:No_Algorithm
    Table_ave(Select_Algorithm,:)=[mean(Table(:,1:5,Select_Algorithm)) std(Table(:,5,Select_Algorithm)) mean(Table(:,6,Select_Algorithm))]
end


%% Baseline and Perfect competition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BL_sol=size(sol_test);
BL_sol=zeros(BL_sol(1),BL_sol(2));
[BL, BL_profits,BL_Penalties, BL_TE] = feval('fitnessFun_WCCI2020',BL_sol,caseStudyData,otherParameters);
%BL
T_BL_profits=sum(BL_profits);
BL_individual_prof=sum(BL_TE.profits,2);

Table_ave(Select_Algorithm+1,:)=[sum(BL_individual_prof) sum(BL_individual_prof(1:25)) sum(BL_individual_prof(26:26)) sum(BL_individual_prof(27:31)) 0 0 0];

%% ___________________________________________________
%Create Figure 3A 
ind_plot=2;
%Figures comparing LM/WS transactions for each agent -9 agent case
[~,Winner_alg]=max(Table_ave(:,1));

figure(ind_plot)
%Joao Suggestion
buy=(-1.*Matrix_energyBuy(:,:,Winner_alg));
sell=(-1.*Matrix_energySell(:,:,Winner_alg));
Balance=buy+sell;
Balance=[sum(Balance(:,1:5),2) sum(Balance(:,6:25),2) sum(Balance(:,27:31),2)] 
bar(Balance','stacked')
 xlabel('Agents (1: Consumers; 2: Prosumers;  3: Generators )')
 ylabel('Energy (kWh)')
 legend('LEM','Agg+WS')
 
 %% Plot clearing prices
WS_CP=WS_price;
Table_CPrices
grid_price=WS_price+Agg_fee;
feed_price=ones(1,length(Table_CPrices(Winner_alg,:)))*Feed_in;
%agg_price=ones(1,T)* Agg_price
load('Results_IJEPS_rev1\CP_bestSol_PT_30AnoPV','CP_bestSol_PT') %It charges Pricetaker solution so PT should be PT=0
%CP_bestSol_PT=Table_CPrices(Winner_alg,:);
X=[grid_price;Table_CPrices(Winner_alg,:);CP_bestSol_PT;feed_price;WS_CP]
CP_indfig=ind_plot+1;
figure(CP_indfig)
plot(X')