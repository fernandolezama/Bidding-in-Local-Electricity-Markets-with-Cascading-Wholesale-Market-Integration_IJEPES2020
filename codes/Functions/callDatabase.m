% available scenarios
function [caseStudyData, DB_name]=callDatabase(scenario)

caseStudyData=[];

switch scenario
    case 1
        load('case_studies\9Agents')
        DB_name='9agents';
        
    case 2
        load('case_studies\5C_20P_5CHP_1day')
        DB_name='5C_20P_5CHP_1day';
         %Line used to increase the PV generation 15 percent more
        caseStudyData.Type2.Gen=caseStudyData.Type2.Gen.*15;
    case 3
        
         load('case_studies\25C_5CHP_1day')
        DB_name='25C_0P_5CHP_1day';
        
    otherwise
        fprintf(1,'Case study not available\n');
        
end


end
