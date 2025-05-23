% run_all_task11.m
patientPaths     = {'1_2','2_7','3_8','4_15'};
consoleEstimates = [9.13, 4.26, 16.29, 10.70];
roiRadius        = sqrt(200/pi);

task_11_LM_model(patientPaths, consoleEstimates, roiRadius);
