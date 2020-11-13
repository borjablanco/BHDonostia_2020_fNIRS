
addpath('/home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/BrainNetViewer_20191031/');


T = importdata('/Users/borjablanco/Downloads/BrainHack_project/BHDonostia_2020_fNIRS/data_files/coords_bn.xlsx');
load /home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/P
U = betas18; %put here the value of each nodes
T(:,5) = U(:,1); 
tmp = char(P);
node = strcat(num2str(T),tmp);
node(:,end-4:end-3) = char(' ');
dlmwrite('/Users/borjablanco/Downloads/BrainHack_project/BHDonostia_2020_fNIRS/data_files/files_for_bnv/plot.node',node,'delimiter','');

cd '/home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/files_for_bnv/'
BrainNet_MapCfg('Akiyama_6Month_3T.nv','plot.node', 'option.mat');

