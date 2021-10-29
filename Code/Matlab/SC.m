% Calculate spatial chaos for all samples
addpath('Spatial Chaos')

%path to folder that containts all the TIL maps placed in different
%folders, please mind that path needs to end with / 
path_to_TIL_folder = '';

% cancer to analyze 
type = 'BRCA';

filePattern = fullfile(append(path_to_TIL_folder,type), '*.png'); 
theFiles = dir(filePattern);
num = length(theFiles);
%path to excel containing information about patients
excel_table = readtable('Drive:\path\TCGA-CDR-SupplementalTableS1.xlsx', 'VariableNamingRule', 'preserve');

% the size is unknown but the algorithm requires a specific format
% we want the data to be vertical, without the 2 it would become horizontal
values = cell(2,1);
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

for k = 1:num
    
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    img = imread(fullFileName);
    img2 = double(img(:,:,1));
    index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
    tabulate_img2 = tabulate(img2(:));
     
    if (isempty(index) == 1)
        fprintf(2, append(baseFileName,' NOT FOUND skiping...\n')); 
        continue
    else
        
        try
            values{counter} = chaos(img2,0);
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

chaos_table_sorted = sortrows(table(names,values,age,race,tumor_stage, ...
    status,os,os_time,DSS,DSS_time, DFI, DSS_time, PFI, PFI_time,patches_num, tll_perc),-2); 

writetable(chaos_table_sorted,append(type,'_chaos_table_sorted.xlsx'))