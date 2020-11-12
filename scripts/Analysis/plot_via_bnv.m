addpath('/home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/BrainNetViewer_20191031/');


T = importdata('/home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/coords_bn.xlsx');
load /home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/P
U = rand(size(P,1),1); %put here the value of each nodes
T(:,5)=U(:,1);
tmp = char(P);
node = strcat(num2str(T),tmp);
node(:,end-4:end-3) = char(' ');
dlmwrite('/home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/plot.node',node,'delimiter','');

cd '/home/davide/Desktop/DAVE/BHDonostia_2020_fNIRS/data_files/files_for_bnv/'
BrainNet_MapCfg('Akiyama_6Month_3T.nv','plot.node', 'option.mat');

