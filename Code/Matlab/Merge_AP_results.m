%% Combine data

% path to excel file generated by R studio
affinity_excel_table = readtable('', 'VariableNamingRule', 'preserve');

%path to excel containing information about patients
excel_table = readtable('Drive:\path\TCGA-CDR-SupplementalTableS1.xlsx', 'VariableNamingRule', 'preserve');
counter = 1;

size = height(affinity_excel_table);
names = cell(size,1);
age= cell(size,1);
race = cell(size,1);
tumor_stage = cell(size,1);
status = cell(size,1);
os = cell(size,1);
os_time = cell(size,1);
DSS = cell(size,1);
DSS_time= cell(size,1);
DFI = cell(size,1);
DFI_time = cell(size,1);
PFI = cell(size,1);
PFI_time = cell(size,1);
patches_num = cell(size,1);
tll_perc = cell(size,1);
c_index = cell(size,1);
ball_hall = cell(size,1);
banfeld_raftery = cell(size,1);
np_mean = cell(size,1);
np_sd = cell(size,1);
wcd_mean = cell(size,1);
wcd_sd = cell(size,1);
cluster_num = cell(size,1);

for k = 1:size
    
     %n = affinity_excel_table.file_names(k);
     n = affinity_excel_table.file_names(k);
     index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(n,1,12));
     
    if (isempty(index) == 1)
        fprintf(2, append(n,' NOT FOUND skiping...\n')); 
        continue
    else
        
        try
            names{counter} = n;
            age{counter} = excel_table.age_at_initial_pathologic_diagnosis(index);
            race{counter} = excel_table.race(index);
            tumor_stage{counter} = excel_table.ajcc_pathologic_tumor_stage(index);
            status{counter} = excel_table.vital_status(index);
            os{counter} = excel_table.OS(index);
            os_time{counter} = excel_table.OS_time(index);
            DSS{counter} = excel_table.DSS(index);
            DSS_time{counter} = excel_table.DSS_time(index);
            DFI{counter} = excel_table.DFI(index);
            DFI_time{counter} = excel_table.DFI_time(index);
            PFI{counter} = excel_table.PFI(index);
            PFI_time{counter} = excel_table.PFI_time(index);
            
            % affinity execl data           
            patches_num{counter} = affinity_excel_table.patches_num(k);
            cluster_num{counter} = affinity_excel_table.clusters_num(k);
            tll_perc{counter} = affinity_excel_table.tll_perc(k);
            c_index{counter} = affinity_excel_table.c_index(k);
            ball_hall{counter} = affinity_excel_table.ball_hall(k);
            banfeld_raftery{counter} = affinity_excel_table.banfeld_raftery(k);
            np_mean{counter} = affinity_excel_table.np_mean(k);
            np_sd{counter} = affinity_excel_table.np_sd(k);
            wcd_mean{counter} = affinity_excel_table.wcd_mean(k);
            wcd_sd{counter} = affinity_excel_table.wcd_mean(k);
                       
            counter = counter + 1;
        catch
            fprintf(2,' FAILED\n');
            continue
        end
    end    
end

affinity_table_sorted = sortrows(table(names,age,race,tumor_stage, ...
    status,os,os_time,DSS,DSS_time, DFI, DSS_time, PFI, PFI_time, patches_num,cluster_num, tll_perc, ...
    c_index, ball_hall, banfeld_raftery, np_mean, np_sd, wcd_mean, wcd_sd),-2); 

writetable(affinity_table_sorted,append(type,'_chaos_table_sorted.xlsx'));