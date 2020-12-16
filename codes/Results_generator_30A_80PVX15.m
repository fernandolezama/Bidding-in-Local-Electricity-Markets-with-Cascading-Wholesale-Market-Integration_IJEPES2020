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
%DB_name='25C_0P_5CHP_1day'
DB_name='5C_20P_5CHP_1day'
PT=0;
DB=2; %1 (500) and 2 (1); %Select the database you want to analyze

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
     algorithm=['Results_IJEPS_rev1\' alg '_' DB_name '_simPV_AGG15WS_PTLM' num2str(PT)]; %'The participants should include their algorithm here'
    
     
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
    
    %Line used to increase the PV generation 15 percent more
    caseStudyData.Type2.Gen=caseStudyData.Type2.Gen.*15;
    
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
        %ResDB(iRuns).Fit_and_p
        %S_val
        profits_a=sum(TE.profits,2);
        M_profitT(iRuns)=sum(profits); %Total profits of the system
        M_fitness(iRuns)=S_val; %Fitness of the system
        M_profitCons(iRuns)=sum(profits_a(1:5,1));
        M_profitProsu(iRuns)=sum(profits_a(6:25,1));
        M_profitProd(iRuns)=sum(profits_a(26:30,1));
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

Table_ave(Select_Algorithm+1,:)=[sum(BL_individual_prof) sum(BL_individual_prof(1:5)) sum(BL_individual_prof(6:25)) sum(BL_individual_prof(26:30)) 0 0 0];

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
Balance=[sum(Balance(:,1:5),2) sum(Balance(:,6:25),2) sum(Balance(:,26:30),2)] 
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
load('Results_IJEPS_rev1\CP_bestSol_PT_30APVx15','CP_bestSol_PT') %It charges Pricetaker solution so PT should be PT=0
%CP_bestSol_PT=Table_CPrices(Winner_alg,:);
X=[grid_price;Table_CPrices(Winner_alg,:);CP_bestSol_PT;feed_price;WS_CP]
CP_indfig=ind_plot+1;
figure(CP_indfig)
plot(X')

% ind_plot=2;
% for Select_Algorithm=1:No_Algorithm %Algorithm distribution of profits by agent
%     figure(ind_plot) %algorithm 1
%     
%     bar(-1.*MH_individual_prof(:,Select_Algorithm))
%     axis([0 10 -1 4])
%     hold on
%     averageRef=repmat(-1.*mean(MH_individual_prof(:,Select_Algorithm)),1,length(MH_individual_prof(:,Select_Algorithm)));
%     plot(averageRef)
%     xlabel('Agents (1 to 3: Consumers; 4 to 6: Prosumers; 7 to 9: Producers)')
%     ylabel('Average cost (EUR)')
%     
%     ind_plot=ind_plot+1;
% end
% 
% %% Plot of the base line
% figure(ind_plot); %Baseline
% bar(-1.*BL_individual_prof)
% axis([0 10 -1 4])
% hold on
% averageRef=repmat(-1.*mean(BL_individual_prof),1,length(MH_individual_prof(:,1)));
% plot(averageRef)
% xlabel('Agents (1 to 3: Consumers; 4 to 6: Prosumers; 7 to 9: Producers)')
% ylabel('Average cost (EUR)')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Figures of the LM and Energy transacted in the WS/retail market
% ind_plot=ind_plot+1;
% for Select_Algorithm=1:No_Algorithm %Algorithm distribution of profits by agent
%     figure(ind_plot)
%     
%     bar(-1.*Matrix_energyBuy(:,1:6,Select_Algorithm)')
%     xlabel('Agents (1 to 3: Consumers; 4 to 6: Prosumers)')
%     ylabel('Energy Bougth(kW)')
%     legend('Local market','Grid')
%         
%     ind_plot=ind_plot+1;
% end
% 
% for Select_Algorithm=1:No_Algorithm %Algorithm distribution of profits by agent
%     figure(ind_plot)
%     
%     bar(-1.*Matrix_energySell(:,4:9,Select_Algorithm)')
%     xlabel('Agents (4 to 6: Prosumers; 7 to 9: Generators )')
%     ylabel('Energy Sell(kWh)')
%     legend('Local market','Grid')
%         
%     ind_plot=ind_plot+1;
% end
% 
%  WS_price=[0.049228	0.047272	0.044936	0.043782	0.043554	0.045226	0.046564	0.048322	0.049872	0.050096	0.049628	0.049796	0.050116	0.050712	0.049648	0.048644	0.048066	0.047712	0.047562	0.048038	0.0491	0.051456	0.050584	0.049478];
%     Original=0.28; %Grid tariff (Can be modified in function of the case study)
%     Feed_in=0.095; %Portugal
%     Agg_fee=0.15; %Considered aggregator fee
%     
% 
% %% Plot clearing prices
% WS_CP=WS_price;
% Table_CPrices
% grid_price=WS_price+Agg_fee;
% feed_price=ones(1,length(Table_CPrices(1,:)))*Feed_in;
% %agg_price=ones(1,T)* Agg_price
% X=[grid_price;Table_CPrices(2,:);feed_price;WS_CP]
% figure(ind_plot)
% plot(X')
% 
% % figure(Select_Algorithm+3) %algorithm 1
% % for i=1:length(profits)
% %     energy_LM(:,i)=TE.TransE(:,1,i); %[LM; GridT; demand/generation]
% %     energy_Grid(:,i)=TE.TransE(:,2,i); %[LM; GridT; demand/generation]
% %     %%** in our LM ISAP paper algorithm negative is buying (demand) and positive is selling (generation)
% %     %%** we change the sign later in line 229
% % end
% % 
% % CON_mask_LM=energy_LM<0; %Buy
% % CON_mask_Grid=energy_Grid<0; %Buy
% % 
% % energy_LM_cum_org=sum(energy_LM.*CON_mask_LM,2);
% % energyG_GEN_org=sum(energy_Grid.*CON_mask_Grid,2);
% % 
% % Matrix_energyPlot=[energy_LM_cum_org';energyG_GEN_org'];
% % bar(-1.*Matrix_energyPlot(:,1:6)')
% % xlabel('Agents (1 to 3: Consumers; 4 to 6: Prosumers)')
% % ylabel('Energy Bougth(kW)')
% % legend('Local market','Grid')
% 
% % No_LM=BL_individual_prof; %profits
% % check=sum(BL_TE.profits,2);
% % LM_only=MH_individual_prof(:,Select_Algorithm)
% 
% 
% % 
% % 
% % for ii=1:3
% %     switch ii
% %         case 1
% %             Success_WS='0.50';
% %         case 2
% %             Success_WS='0.60'; %.75 is not giving the results expecting accepting almost all bids always
% %         case 3
% %             Success_WS='1';
% %     end
% %     file_name=[algorithm '_' Success_WS 'WS.mat'];
% %     load(file_name)%, 'WS_Results', 'WS_profits', 'Total_WS', 'time_step','Consumer_cost_cum',)
% %     
% %     Agg_price=0.2;
% %     
% %     Success_mask=CP-0.05==0
% %     
% %     Agg_cost=(-1*(energyG_matrix.*CON_mask)).*repmat(CP-0.05,length(energyG_matrix(:,1)),1);
% %     Agg_cost_cum=sum(Agg_cost,2)
% %     
% %     Current_cost=(-1*(energyG_matrix.*CON_mask)).*repmat(ones(1,T)*0.28,length(energyG_matrix(:,1)),1);
% %     Current_cost_cum=sum(Current_cost,2)
% %     
% %     Consumer_cost=(-1*(energyG_matrix.*CON_mask)).*repmat((~Success_mask*Agg_price)+(Success_mask*0.28),length(energyG_matrix(:,1)),1);
% %     Consumer_cost_cum=sum(Consumer_cost,2)
% %     
% %     Consumer_cost_onlyWS=(-1*(energyG_matrix.*CON_mask)).*repmat((~Success_mask*Agg_price),length(energyG_matrix(:,1)),1);
% %     Consumer_cost_onlyWS_cum=sum(Consumer_cost_onlyWS,2)
% %     
% %     Consumer_cost_restWS=(-1*(energyG_matrix.*CON_mask)).*repmat((Success_mask*0.28),length(energyG_matrix(:,1)),1);
% %     Consumer_cost_restWS_cum=sum(Consumer_cost_restWS,2)
% %     
% %     
% %     Consumer_reduction=Current_cost-Consumer_cost;
% %     Aggregator_profit=Consumer_cost_onlyWS-Agg_cost;
% %     
% %     LM_WS(:,ii)=LM_only-sum(-1.*Consumer_reduction,2)
% %     
% %     figure(Select_Algorithm+3+ii) %algorithm 1
% %     energy_WS_cum=energyG_matrix.*CON_mask.*repmat(~Success_mask,length(energyG_matrix(:,1)),1)
% %     
% %     Matrix_energyPlot=[energy_LM_cum_org';energyG_GEN_org'-sum(energy_WS_cum,2)';sum(energy_WS_cum,2)'];
% %     
% %     
% %     
% %     bar1=bar(-1.*Matrix_energyPlot(:,1:6)','stacked','FaceColor','flat')
% %     set(bar1(3),'DisplayName','Local market',...
% %         'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461]);
% %     set(bar1(2),'DisplayName','Grid',...
% %         'FaceColor',[0.501960813999176 0.501960813999176 0.501960813999176]);
% %     set(bar1(1),'DisplayName','Wolesale',...
% %         'FaceColor',[0.24705882370472 0.24705882370472 0.24705882370472]);
% %     %     b.CData(1,:) = [0.831372559070587 0.815686285495758 0.7843137383461];
% %     %      b.CData(2,:) = [0.501960813999176 0.501960813999176 0.501960813999176];
% %     %       b.CData(3,:) = [0.24705882370472 0.24705882370472 0.24705882370472];
% %     xlabel('Agents (1 to 3: Consumers; 4 to 6: Prosumers)')
% %     ylabel('Energy (kWh)')
% %     legend('Local market','Grid','Wholesale')
% %     
% %     Aggregator(ii,:)=[sum(sum(Aggregator_profit)) sum(Agg_cost_cum) sum( Consumer_cost_restWS_cum) sum( Consumer_cost_onlyWS_cum) sum( Current_cost_cum)] ;
% %     
% % end
% % 
% % Table_1=[No_LM' sum(No_LM); LM_only' sum(LM_only); [LM_WS' sum(LM_WS)']]
% % 
% % figure(Select_Algorithm+7) %algorithm 1
% % 
% % bar(Aggregator(:,1:3),'stacked')
% % xlabel('Different SR in the WS')
% % ylabel('Cost (m.u)')
% % legend('Aggregator profit','Aggregator WS cost')
% % hold on
% % plot(Aggregator(:,5))
% % legend('Consumer cost without Aggregator')
% % 




