%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:          Fernando Lezama
% Modified by FLC \GECAD 2019

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

function [Fit_and_p,FVr_bestmemit, fitMaxVector, bestFitness ] = ...
    ACO_LM(acoParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit)

S_bestval=Inf;
FVr_bestmemit=0;

%% Initial parameters of ACO
maxIter = acoParameters.maxIter; %500
antNo = acoParameters.antNo; %50
rho = acoParameters.rho; %0.5; % Evaporation rate
alpha = acoParameters.alpha; %1;  % Phromone exponential parameters
beta = acoParameters.beta; %1;  % Desirability exponetial paramter

caseStudyData.low_habitat_limit=low_habitat_limit;
I_D          = numel(up_habitat_limit); %Number of variables or dimension
acoParameters.nVariables=I_D;
FVr_minbound = low_habitat_limit;
FVr_maxbound = up_habitat_limit;
%I_itermax    = acoParameters.maxIter; %500
fnc= otherParameters.WCCI_2020_funct;

numType1=caseStudyData.General.numN1; %Consumers only
numType2=caseStudyData.General.numN2; %Prosumers
%numType3=caseStudyData.General.numN3; %Consumers with Battery
numType4=caseStudyData.General.numN4; %Generators
N_A=numType1+numType2+numType4;
T=caseStudyData.General.numPeriods;
DiscretLevel=acoParameters.DiscretLevel;

%% tau initialization
tau0 = 1/DiscretLevel; %10 * 1 / (  graph.n * mean( graph.edges(:)  )  );  % Initial phromone concentration
tau.bid = tau0 * ones( T , DiscretLevel, N_A); % Phromone matirx
tau.price = tau0 * ones( T , DiscretLevel, N_A); % Phromone matirx

%% Calculate bidding discretization
bid_lb=FVr_minbound(1,1:I_D/2);
price_lb=FVr_minbound(1,I_D/2+1:end);
bid_ub=FVr_maxbound(1,1:I_D/2);
price_ub=FVr_maxbound(1,I_D/2+1:end);
A_bid_lb=reshape(bid_lb,T,N_A);
A_price_lb=reshape(price_lb,T, N_A);
A_bid_ub=reshape(bid_ub,T,N_A);
A_price_ub=reshape(price_ub,T, N_A);
Edges.bid=zeros( T , DiscretLevel, N_A); %Possible actions to choose by agents
Edges.price=zeros( T , DiscretLevel, N_A); %Possible actions to choose by agents
for i=1:N_A
    for j=1:T
        %Quantity
        lb_Q=A_bid_lb(j,i);
        ub_Q=A_bid_ub(j,i);
        %step=linspace(lb_Q, ub_Q, DiscretLevel)
        Edges.bid(j,:,i)=linspace(lb_Q, ub_Q, DiscretLevel);
        
        %price
        lb_P=A_price_lb(j,i);
        ub_P=A_price_ub(j,i);
        %step=linspace(lb_Q, ub_Q, DiscretLevel)
        Edges.price(j,:,i)=linspace(lb_P, ub_P, DiscretLevel);
    end
end

%% eta = 1./ ones( T , DiscretLevel, N_A);  % desirability of each edge
% %% This was the form used for EEM paper
% rhoma=0.1*(1:1:DiscretLevel); %This is a factor for the pheromone of desirability
% for i=1:N_A
%     if i<=numType1
%         eta.bid(:,:,i)=Edges.bid(:,:,i)+rhoma;
%         eta.price(:,:,i)=1./Edges.price(:,:,i); %We use the inverse because consumers want to buy to a smaller price
%     end
%     if i>numType1 && i<=numType1+numType2
%         gen_ind=sum(Edges.bid(:,:,i)<0,2)>0;
%         eta.bid(~gen_ind,:,i)=Edges.bid(~gen_ind,:,i)+rhoma; %Positive values try to buy as much as possible
%         eta.bid(gen_ind,:,i)=abs(Edges.bid(gen_ind,:,i)+flip(rhoma)); %Negative values try to sell as much as possible
%         %eta.bid(:,:,i)=abs(Edges.bid(:,:,i)+rhoma); %Negative values to positive because they prefer to sell as much as possible
%         eta.price(~gen_ind,:,i)=1./Edges.price(~gen_ind,:,i);
%         eta.price(gen_ind,:,i)=Edges.price(gen_ind,:,i)+rhoma;
%     end
%     if i>numType1+numType2 && i<=N_A
%         eta.bid(:,:,i)=abs(Edges.bid(:,:,i))+flip(rhoma);
%         eta.price(:,:,i)=Edges.price(:,:,i)+rhoma; %We desire to sell to the higest price price
%     end
% end
%% This is an alternative way of creating the eta matrix

Vec_eta=cumsum(ones(1,DiscretLevel)./DiscretLevel);
%Vec_eta=ones(1,DiscretLevel);

favors_ascendt=1./repmat(flip(Vec_eta),T,1);
favors_descendt=1./repmat(Vec_eta,T,1); 
for i=1:N_A
    if i<=numType1
        eta.bid(:,:,i)=favors_ascendt; %Consumers favor to buy in the market
        eta.price(:,:,i)=favors_descendt; %We use the inverse because consumers want to buy to a smaller price
    end
    if i>numType1 && i<=numType1+numType2
        gen_ind=sum(Edges.bid(:,:,i)<0,2)>0;
        
        eta.bid(~gen_ind,:,i)=favors_ascendt(~gen_ind,:); %Positive values try to buy as much as possible
        eta.price(~gen_ind,:,i)=favors_descendt(~gen_ind,:); %Consumers want to buy at cheapest price
        
        eta.bid(gen_ind,:,i)=favors_descendt(gen_ind,:); %Negative values try to sell as much as possible
        eta.price(gen_ind,:,i)=favors_ascendt(gen_ind,:);
    end
    if i>numType1+numType2 && i<=N_A
        eta.bid(:,:,i)=favors_descendt;
        eta.price(:,:,i)=favors_ascendt; %We desire to sell to the higest price price
    end
end

graph.n=N_A;
graph.Edges=Edges;

%% Main loop of ACO
bestFitness = inf*ones(1,N_A);
bestTour = ones(1,N_A); %Index of best ant from solutions

for t = 1 : maxIter
    % Create Ants
    colony = [];
    [colony, convergence_flag]= createColony(graph, colony , antNo, tau, eta, alpha,  beta);
    %[colony, convergence_flag] = createColony_dym(graph, colony , antNo, tau, eta, alpha,  beta, caseStudyData);
    if any(convergence_flag(:) == 1)
        [x1,y1,z1] = ind2sub(size(convergence_flag),find(convergence_flag == 1));
        for i_corr=1:length(x1)
            if y1==1
                 tau.bid(x1(i_corr),:,z1(i_corr))=tau0 ;
            else
                 tau.price(x1(i_corr),:,z1(i_corr))=tau0;
            end
        end
        %% Drastic restart pheromone matrix
%          for i_corr=1:length(x1)
%             if y1==1
%                  tau.bid(x1(i_corr),:,:)=tau0 ;
%             else
%                  tau.price(x1(i_corr),:,:)=tau0;
%             end
%         end
    end
    
    
    %% Decoding ant information in the form of MH solution
    for i=1:antNo
        [~,N_A]=size(colony);
        sol_q=[];
        sol_p=[];
        for j=1:N_A
            ind_bid=colony(j).ant(i).bid; %Index of solutions
            ind_price=colony(j).ant(i).price; %Index of solutions
           A=graph.Edges.bid(:,:,j);
            B=graph.Edges.price(:,:,j);
            for tt=1:T
                q(tt)=A(tt, ind_bid(tt));
                p(tt)=B(tt, ind_price(tt));
            end
            sol_q=[sol_q  q];
            sol_p=[sol_p  p];
        end
        Solutions(i,:)=[sol_q,sol_p];
    end
    % Calculate the fitness of all the colony (across agents and across ants
    [S_val, Other,penalties,TE]=feval(fnc, Solutions,caseStudyData,otherParameters);
    
    %% Store best solution found
    [S_bestval_temp,I_best_index] = min(S_val);
    if S_bestval_temp<S_bestval
        FVr_bestmemit = Solutions(I_best_index,:); % best member of current iteration
        fitMaxVector(1,t) = S_bestval_temp;
        S_bestval=S_bestval_temp;
    else
        FVr_bestmemit = FVr_bestmemit;
        fitMaxVector(1,t)=S_bestval;
    end
    
    
    
    %% Calculate the fitness values of all ants individually
    for j=1:N_A
        for i = 1 : antNo
            colony(j).ant(i).fitness = -1.*TE(i).profits(j,:)'; %Already in the form of cost for minimization
            %% Constraints Include pelaties in the learning process
            colony(j).ant(i).P_agent=penalties(j,:,i)';
        end
    end
    
    %% Find the best ant (queen) %Best local for each agent
    for j=1:N_A
        %for i = 1 : antNo
        %allAntsFitness(j,:)=sum([colony(j).ant(:).fitness]);
        %% Constraints
        allAntsFitness(j,:)=sum([colony(j).ant(:).fitness])+sum([colony(j).ant(:).P_agent]);
        
        [ minVal , minIndex ] = min( allAntsFitness(j,:) );
        
        if minVal < bestFitness(1,j)
            bestFitness(t,j) = minVal; %   colony.ant(minIndex).fitness;
            bestTour(t,j) = minIndex;
            Tfitnes(j,:,t)=colony(j).ant(minIndex).fitness;
            Tpenal(j,:,t)=colony(j).ant(minIndex).P_agent;
        else
            bestFitness(t,j) = bestFitness(t-1,j); %   colony.ant(minIndex).fitness;
            bestTour(t,j) = bestTour(t-1,j);
            Tfitnes(j,:,t)=Tfitnes(j,:,t-1);
            Tpenal(j,:,t)=Tpenal(j,:,t-1);
        end
        colony(j).queen.tour = bestTour(t,j);
        colony(j).queen.bid=colony(j).ant(bestTour(t,j)).bid;
        colony(j).queen.price=colony(j).ant(bestTour(t,j)).price;
        colony(j).queen.bestFitness = bestFitness(t,j);
        colony(j).queen.Tfitness(:,1)=Tfitnes(j,:,t)';
        colony(j).queen.Tpenal(:,1)=Tpenal(j,:,t)';
    end
   
    %   for j=1:N_A
    %     colony(j).queen.tour = bestTour(t,j);
    %     colony(j).queen.bid=colony(j).ant(bestTour(t,j)).bid;
    %     colony(j).queen.price=colony(j).ant(bestTour(t,j)).price;
    %     colony(j).queen.fitness = bestFitness(t,j);
    %
    %     %Decoding best so far solution
    %     ind_bid=colony(j).queen.bid; %Index of solutions
    %     ind_price=colony(j).queen.price; %Index of solutions
    %
    %     A=graph.Edges.bid(:,:,j);
    %     B=graph.Edges.price(:,:,j);
    %            for tt=1:T
    %               q(tt)=A(tt, ind_bid(tt));
    %               p(tt)=B(tt, ind_price(tt));
    %
    %            end
    %           sol_q=[sol_q  q];
    %           sol_p=[sol_p  p];
    %   end
    %    FVr_bestmemit = [sol_q,sol_p];
    %
    
    %% Update phromone matrix
    tau = updatePhromone(tau , colony, caseStudyData );
    
    %% Evaporation
    %heatmap(tau.bid(:,:,8))
    %pause(1)    
    tau.bid= ( 1 - rho ) .* tau.bid;
    tau.price= ( 1 - rho ) .* tau.price;
    % heatmap(tau.bid(:,:,8))
    %pause(1)
   
    %     % Display the results
    %     %% update best results
    %     [S_bestval,I_best_index] = min(S_val);
    %     %FVr_bestmemit = Solutions(I_best_index,:); % best member of current iteration
    %      % store fitness evolution and obj fun evolution as well
    %     fitMaxVector(1,t) = S_bestval;
    %S_bestval
    fprintf('Fitness value: %f\n',S_bestval )
    fprintf('Generation: %d\n',t)
    
end

Fit_and_p=S_bestval;
% figure (1), plot(bestFitness(:,1:3))
% figure (2), plot(bestFitness(:,4:6))
% figure (3), plot(bestFitness(:,7:9))



function [ tau ] = updatePhromone(tau , colony, caseStudyData)
Queen_flag=0;

% Update the phromone matrix.
%numType1=caseStudyData.General.numN1; %Consumers only
%numType2=caseStudyData.General.numN2; %Prosumers
[~,N_A]=size( colony);
[~,antNo]=size( colony(1).ant);
%eps=1e-6;
%rhoma=0.1;
%for j=1:N_A
for j=1:N_A %Only reinforce Generators behaivour
    for i = 1 : antNo % for each ant
       % ant_cont=zeros(length(colony(j).ant(i).fitness(:)),1);
        %FIT=colony(j).ant(i).fitness(:);
        %% Constraint
        FIT=colony(j).ant(i).fitness(:)+colony(j).ant(i).P_agent(:);
%         ind_gen=(FIT)<0;
%         ind_con=(FIT)>0;
%         ind_neutro=(FIT)==0;
        %% Why I used an inverse function for this
        %x1=-5:0.5:5
        %x2=0.01:0.05:10
        %y1=1./(exp(1).^x1)
        %y2=1./x2
        %plot(x1,y1)
        %hold on
        %plot(x2,y2)
        
        
        %% Since the fitness can take negative values, we have three cases %Inverse of the Fitness (More pheromone to better solutions
%         ant_cont( ind_con)=1./(FIT(ind_con)+eps); %C1 Inverse when fitness is positive (reinforce the minimum)
%         ant_cont( ind_gen)=abs(FIT(ind_gen)+eps)*rhoma; %C2 abs multiplied by a factor (reinforce the minimum negative)
%         ant_cont(ind_neutro)=0; % C3 if 0, do not reinforce any solution
%         if j<=numType1
%             ant_cont=1./(FIT+eps);
%         else
%             ant_cont=1./(exp(1).^FIT);
%         end
           ant_cont=1./(exp(1).^FIT);
           %ant_cont=1./(FIT);
           %ant_cont=1./100.^FIT;
         
        ant_bid_ind=colony(j).ant(i).bid;
        ant_price_ind=colony(j).ant(i).price;
        
        for t=1:length(ant_bid_ind)
            tau.bid(t,ant_bid_ind(t),j)=tau.bid(t,ant_bid_ind(t),j)+ant_cont(t); %Update of bid pheromone
            tau.price(t,ant_price_ind(t),j)=tau.price(t,ant_price_ind(t),j)+ant_cont(t); %Update price pheromone
        end
    end
    
    %% Adding extrapheromone for queen Ant
    %Qant_cont=zeros(length(colony(j).ant(i).fitness(:)),1);
    %FIT=colony(j).ant(i).fitness(:);
    %% Constraint
    if Queen_flag==1
        QFIT=colony(j).queen.Tfitness+ colony(j).queen.Tpenal;
        Qant_cont=1./(exp(1).^QFIT);

        Qant_bid_ind=colony(j).queen.bid;
        Qant_price_ind=colony(j).queen.price;
        for t=1:length(Qant_bid_ind)
            tau.bid(t,Qant_bid_ind(t),j)=tau.bid(t,Qant_bid_ind(t),j)+Qant_cont(t); %Update of bid pheromone
            tau.price(t,Qant_price_ind(t),j)=tau.price(t,Qant_price_ind(t),j)+Qant_cont(t); %Update price pheromone
        end
    end
end


function [ colony, tau_reset_flag ] = createColony( graph, colony , antNo, tau, eta, alpha,  beta)
[~,~,N_A]=size(graph.Edges.bid); %Number of agents

for a=1:N_A
    tau_temp_bid=tau.bid(:,:,a);
    tau_temp_price=tau.price(:,:,a);
    eta_temp_bid=eta.bid(:,:,a);
    eta_temp_price=eta.price(:,:,a);
    
     b_P_allNodes = tau_temp_bid .^ alpha .* eta_temp_bid.^ beta;
     b_P = b_P_allNodes ./ sum(b_P_allNodes,2);
    
     p_P_allNodes = tau_temp_price .^ alpha .* eta_temp_price.^ beta;
     p_P = p_P_allNodes ./ sum(p_P_allNodes,2);
     
    for i = 1 : antNo
        nextNode = rouletteWheel(b_P);
        colony(a).ant(i).bid=nextNode;
        
        nextNode = rouletteWheel(p_P);
        colony(a).ant(i).price=nextNode;
    end
   tau_reset_flag(:,:,a)=[sum(b_P==1,2),sum(p_P==1,2)];
end


function [ colony, tau_reset_flag ] = createColony_dym( graph, colony , antNo, tau, eta, alpha,  beta,caseStudyData)

%% This part is needed for correction
[T,Level,N_A]=size(graph.Edges.bid); %Number of agents
numType1=caseStudyData.General.numN1; %Consumers only
numType2=caseStudyData.General.numN2; %Prosumers
ind_price=numel(caseStudyData.low_habitat_limit)/2;
% Q_lead=-1*FM_pop(:,T*(numType1+numType2)+1:ind_price);
b_chp=caseStudyData.Type4.MC; % cost factor for the CHP generation
c_CHP=@(b_chp,g_chp) ((b_chp.*sqrt(g_chp))./g_chp+1e-10);

for a=1:N_A
    tau_temp_bid=tau.bid(:,:,a);
    tau_temp_price=tau.price(:,:,a);
    eta_temp_bid=eta.bid(:,:,a);
    eta_temp_price=eta.price(:,:,a);
    
    p_P_allNodes = tau_temp_price .^ alpha .* eta_temp_price.^ beta;
    p_P = p_P_allNodes ./ sum(p_P_allNodes,2);
    b_P_allNodes = tau_temp_bid .^ alpha .* eta_temp_bid.^ beta;
    b_P = b_P_allNodes ./ sum(b_P_allNodes,2);
    %% %%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : antNo
        nextNode = rouletteWheel(b_P);
        colony(a).ant(i).bid=nextNode;
        
        if a>numType1+numType2
            A=graph.Edges.bid(:,:,a);
            B=graph.Edges.price(:,:,a);
            %% This part is needed for correction
            for tt=1:T
                q(tt)=A(tt,  nextNode(tt));
            end
            Q_lead=-1.*q;
            Marginal_Cost=c_CHP(b_chp(1,1),Q_lead);
            minPossibleBid=caseStudyData.low_habitat_limit(1,ind_price+1+T*(a-1):ind_price+T*(a));
            Allow_minbid=max(Marginal_Cost,minPossibleBid);
            
            %inhibiting pheromone
            ind_allow=B> repmat(Allow_minbid',1,Level);
            p_P_allNodes= ind_allow.*p_P_allNodes;
            
            p_P = p_P_allNodes ./ sum(p_P_allNodes,2);
            ind_to_zero=isnan(p_P);
            %correcting 0s
            p_P(ind_to_zero)=0;
            nextNode = rouletteWheel(p_P);
            colony(a).ant(i).price=nextNode;
            
            %% Do i really need to correct the bid?? yes you need to make 0 the energy sell if the marginal cost is too high
            colony(a).ant(i).bid(sum(ind_to_zero,2)>0)=Level; %This is tu put the bid to 0 energy when the marginal cost is higher than the max
        else
            nextNode = rouletteWheel(p_P);
            colony(a).ant(i).price=nextNode;
        end
    end
    
    tau_reset_flag(:,:,a)=[sum(b_P==1,2),sum(p_P==1,2)];
end


function [ nextNode ] = rouletteWheel( P )
% Roulette wheel to choose one edge based on P values
cumsumP = cumsum(P,2);
r = rand(length(P(:,1)),1);
[~, nextNode] = max(r <= cumsumP,[],2);
%nextNode = find( r <= cumsumP );
%nextNode = nextNode(1);
return

