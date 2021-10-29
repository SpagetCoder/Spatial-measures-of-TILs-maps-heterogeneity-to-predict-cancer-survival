%% GLCM on full dataset

%path to folder that containts all the TIL maps placed in different folders
path_to_TIL_folder = '';

% cancer type to analyze 
type = 'BRCA';

filePattern = fullfile(append(path_to_TIL_folder,type), '*.png'); 
theFiles = dir(filePattern);
num2 = length(theFiles);
% path to the excel file that containts information about patiets vital
% status etc.
excel_table = readtable('Drive:\path\TCGA-CDR-SupplementalTableS1.xlsx', 'VariableNamingRule', 'preserve');

% the size is unknown but the algorithm requires a specific format
% we want the data to be vertical, without the 2 it would become horizontal
struct_ness = cell(2,1);
names = cell(2,1);
age = cell(2,1);
race = cell(2,1);
tumor_stage = cell(2,1);
status = cell(2,1);
os = cell(2,1);
os_time = cell(2,1);
DSS = cell(2,1);
DSS_time = cell(2,1);
DFI = cell(2,1);
DFI_time = cell(2,1);
PFI = cell(2,1);
PFI_time = cell(2,1);
patches_num = cell(2,1);
tll_perc = cell(2,1);

counter = 1;

for i = 1:num2
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(theFiles(i).folder, baseFileName);
    img = imread(fullFileName);
    img2 = double(img(:,:,1));
    index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
    tabulate_img2 = tabulate(img2(:));
    
    if (isempty(index) == 1)
        bigboystr = append(baseFileName,' NOT FOUND skiping...');
        disp(bigboystr);
        continue
    else
        
        try
            
        gcm = graycomatrix(img2,'NumLevels',2,'Offset',[0 1; -1 1; -1 0; -1 -1],'Symmetric',false);
        gcm1 = gcm(:,:,1);
        gcm2 = gcm(:,:,2);
        gcm3 = gcm(:,:,3);
        gcm4 = gcm(:,:,4);
        
        struct_ness{counter,1} = gcm1(1,1) + gcm2(1,1) + gcm3(1,1) + gcm4(1,1);
        struct_ness{counter,2} = gcm1(2,2) + gcm2(2,2) + gcm3(2,2) + gcm4(2,2);
        
        names{counter} = baseFileName;
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
        tll_perc{counter} = tabulate_img2(2,3);
        patches_num{counter} = tabulate_img2(2,2);
          
        counter = counter + 1;
        
        catch
            fprintf(2,append(baseFileName,' FAILED TO CALCULATE SC skiping...\n'));
            continue
        end
    end  
end

%% Create table
tab = horzcat(struct_ness);
tab(:,1) = num2cell(cell2mat(tab(:,1))/max(cell2mat(tab(:,1))));
tab(:,2) = num2cell(cell2mat(tab(:,2))/max(cell2mat(tab(:,2))));
tab(:,1) = num2cell(cell2mat(tab(:,1))+cell2mat(tab(:,2)));
co_occurance_sorted_table = splitvars(sortrows(table(names, tab, age, race, tumor_stage, ...
status,os, os_time, DSS, DSS_time, DFI, DFI_time, PFI, PFI_time, ...
patches_num, tll_perc), -2), 'tab', 'NewVariableNames', {'struct_ness_1', 'struct_ness_2'});

writetable(co_occurance_sorted_table,append(type,'_co_occurance_table_sorted.xlsx'))