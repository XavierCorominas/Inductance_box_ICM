


close all 
clear all 
addpath 'D:\DOCUMENTS\LI_rTMS_project\Models_final\coil design\Personal_coils\Manon_coil'

path_to_digimouse =(['D:\DOCUMENTS\LI_rTMS_project\Models_final\coil design\Personal_coils\simulation tests\'])
addpath path_to_digimouse

addpath ('C:\Users\XAVIER\SimNIBS-4.0')

path_to_coil = (['D:\DOCUMENTS\LI_rTMS_project\Models_final\coil design\Personal_coils\Manon_coil']);
addpath path_to_coil

%
% create cells representing coil coordinates and correct for errors
T = readtable('coilCoordinates1mm.xlsx');
H= table2array(T);

G = H;
idf = isnan(G);
find(idf)
G(idf) = [];
G = reshape(G,[],3)

%% 
cd (['D:\DOCUMENTS\LI_rTMS_project\Models_final\coil design\Personal_coils\Manon_coil'])
coil = 'coil_manon.xlsx';
writematrix(G,coil,'Sheet',1,'Range','D1')


%% Check field distribution (visualize magnetic field vector distribution)
app.hdr = nifti_load('nifti_blank.nii');
            
            app.hdr.sform(3,4)=-85;

            [app.xx,app.yy,app.zz]=meshgrid(linspace(-((app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4)),(app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4),app.hdr.dim(2)),...
                                            linspace(-((app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4)),(app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4),app.hdr.dim(3)),...
                                            linspace((app.hdr.dim(4)-1)*app.hdr.sform(3,3)+app.hdr.sform(3,4),app.hdr.sform(3,4),app.hdr.dim(4)));


hdr = nifti_load('coil_manon.nii')

                   n=5; %set resolution
                
q=quiver3(app.xx(1:n:end,1:n:end,1:n:end),app.yy(1:n:end,1:n:end,1:n:end),app.zz(1:n:end,1:n:end,1:n:end),hdr.vol(1:n:end,1:n:end,1:n:end,2),hdr.vol(1:n:end,1:n:end,1:n:end,1),hdr.vol(1:n:end,1:n:end,1:n:end,3));
        

%% Head reconstruction of the rat model 
subject_id = 'rat'
path_to_mri_nifti = ['C:\Users\XAVIER\OneDrive\Escritorio\manon_test\'];
filename = (['WHS_SD_rat_T2star_v1.01.nii.gz']);

cd(path_to_mri_nifti)
clc;
[status,cmdout] = system(['cd ',path_to_mri_nifti],'-echo');

[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\charm',' ',subject_id, ' ' ,filename, ' --forceqform --forcerun'], '-echo');
%[status,cmdout] = system(['/Users/tuomasmutanen/Applications/SimNIBS-4.0/bin/charm',' ',subject_id, ' ' ,filename], '-echo');

%%  Run a simple simulation 
save_folder = [path_to_digimouse,'simulation\']
if ~exist(save_folder)
      mkdir([save_folder]);
end



% COMPUTE OPTIMIZED COIL POSITION
tms_opt = opt_struct('TMSoptimize');
% Select the head mesh
tms_opt.fnamehead = ([path_to_digimouse,'\Mouse_Digimouse.msh']);
% Select output folder
tms_opt.pathfem = (path_to_digimouse);
% Select the coil model
tms_opt.fnamecoil = ([path_to_coil,'\coil_manon.nii']);
% Select a target for the optimization
tms_opt.target = ([-3.5, 4.74, 29.60]);
tms_opt.distance = 2;
tms_opt.target_size = 3;
tms_opt.didt = 1.0e6;
tms_opt.search_radius = 20;
opt.target_direction = [];
tms_opt.solver_options = 'pardiso';
run_simnibs(tms_opt);




% RUN SIMULATION
%Initialize  simulation session
s = sim_struct('SESSION');
s.fnamehead = ([path_to_digimouse,'\Mouse_Digimouse.msh']);
s.fields = 'eEjJ'
% Output folder
s.pathfem = ([path_to_digimouse,'\sim2']);
s.poslist{1} = sim_struct('TMSLIST');
% Select coil
s.poslist{1}.fnamecoil = ([path_to_coil,'\coil_manon.nii']);
%  1 mm distance between coil and head
% slect coil position 
s.poslist{1}.pos;
s.poslist{1}.pos.matsimnibs = [
  -0.9956 -0.0868 -0.035   1.4866
  0.      0.3742 -0.9274  7.4164
  0.0936 -0.9233 -0.3725 30.9926
  0.      0.      0.      1.    ];
s.poslist{1}.pos.distance = 2;
s.poslist{1}.pos.didt = 1.0e6;  % select desired intensity according to you dI/dt values in A/s
% Select coil centre

run_simnibs(s)




%% visualize electric field vectors of different coils as example 





Manon = nifti_load('coil_manon.nii');

fig_butterfly = nifti_load('Magstim_70mm_Fig8.nii.gz');


%% 
app.hdr = nifti_load('nifti_blank.nii');
            
 app.hdr.sform(3,4)=-85;

            [app.xx,app.yy,app.zz]=meshgrid(linspace(-((app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4)),(app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4),app.hdr.dim(2)),...
                                            linspace(-((app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4)),(app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4),app.hdr.dim(3)),...
                                            linspace((app.hdr.dim(4)-1)*app.hdr.sform(3,3)+app.hdr.sform(3,4),app.hdr.sform(3,4),app.hdr.dim(4)));


hdr = nifti_load('coil_manon.nii')

                   n=3; %set resolution
figure(1)
q=quiver3(app.xx(1:n:end,1:n:end,1:n:end),app.yy(1:n:end,1:n:end,1:n:end),app.zz(1:n:end,1:n:end,1:n:end),hdr.vol(1:n:end,1:n:end,1:n:end,2),hdr.vol(1:n:end,1:n:end,1:n:end,1),hdr.vol(1:n:end,1:n:end,1:n:end,3));



app2 = app;

app2.xx(122:end,:,:) = []
app2.xx(:,82:end,:) = []
app2.xx(:,:,42:end) = []

app2.yy(122:end,:,:) = []
app2.yy(:,82:end,:) = []
app2.yy(:,:,42:end) = []

app2.zz(122:end,:,:) = []
app2.zz(:,82:end,:) = []
app2.zz(:,:,42:end) = []



figure(2)
q=quiver3(app2.xx(1:n:end,1:n:end,1:n:end),app2.yy(1:n:end,1:n:end,1:n:end),app2.zz(1:n:end,1:n:end,1:n:end),fig_butterfly.vol(1:n:end,1:n:end,1:n:end,2),fig_butterfly.vol(1:n:end,1:n:end,1:n:end,1),fig_butterfly.vol(1:n:end,1:n:end,1:n:end,3));




%% 

T = table2cell(T)

C(:,1) = T(:,1)

for i = 1:length(T)
    C{i,2} = T{i,2}(end-1:end);
end 


for i = 1:length(T)
    C{i,3} = T{i,3}(end-1:end);
end 

% delate ; charachters
for i = 1:length(C)
  C{i,2} = regexp(C{i,2},';', 'split');
  %emptyCells = cellfun('isempty', {i,2});
 % C{i,2} = []
end 

for i = 1:length(C)
    C{i,2} = join(C{i,2});
 end 


 for i = 1:length(C)
  C{i,3} = regexp(C{i,3},';', 'split');
  %emptyCells = cellfun('isempty', {i,2});
 % C{i,2} = []
end 
for i = 1:length(C)
    C{i,3} = join(C{i,3});
 end 


%%

D= C(:,1);
D = cell2mat(D);


V = C(:,2);
for i = 1:length(V)
    V{i,1} = cell2mat(V{i,1});
end 
%V= strrep(V,' ','');
V = cell2mat(V);
V = str2num(V);



K = C(:,3);
for i = 1:length(K)
    K{i,1} = cell2mat(K{i,1});
end 
K=strrep(K,' ','')
K = cell2mat(K);
K = str2num(K);

% smooth eleiminating outliwers
D = smoothdata(D,'rlowess');
V = smoothdata(V,'rlowess');
K = smoothdata(K,'rlowess');


%coil_manon = 0;
coil_manon(:,1) = D;
coil_manon(:,2) = V;
coil_manon(:,3) = K;

figure(1)
plot3(coil_manon(:,1),coil_manon(:,2),coil_manon(:,3),'-o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')


figure(2)
plot3(coil_manon(:,1),coil_manon(:,2),coil_manon(:,3),'o');





%%
cd 'D:\DOCUMENTS\LI_rTMS_project\Models_final\coil design\Personal_coils\Manon_coil'
T = array2table(coil_manon)
coil = 'coil_manon.xlsx';
writematrix(H,coil,'Sheet',1,'Range','D1')

save coil_manon.mat

%% 
%OTHER TESTS FOR OTHER COILS --> first load the coil

figure(3)
plot3(mycoil(:,1),mycoil(:,2),mycoil(:,3),'-o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')

figure(4)
plot3(mycoil(:,1),mycoil(:,2),mycoil(:,3),'o');

plot3(mycoil_braided_opt(:,4),mycoil_braided_opt(:,5),mycoil_braided_opt(:,6),'-o');

plot3(H(:,1),H(:,2),H(:,3),'-o');



%% 
% EXAMPLE CONTINUOUS SINUSOIDAL TESTS

t = 0:pi/50:10*pi;
st = sin(t)
ct = cos (t)



plot3(st,ct,t,'-o');




%% 






close all 
clear all 
addpath 'D:\DOCUMENTS\LI_rTMS_project\Models_final\coil design\Personal_coils\Manon_coil'

%
% create cells representing coil coordinates
T = readtable('coilCoordinates1.xlsx');
F= table2array(T);
Pcenter = H;
coil_manon = Pcenter;

figure(1)
plot3(coil_manon(:,1),coil_manon(:,2),coil_manon(:,3),'-o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')


figure(2)
plot3(coil_manon(:,1),coil_manon(:,2),coil_manon(:,3),'o');

%%
I0 = 1000;  
dIdt = 2*pi*1e4*I0; 
%dIdt = 0; %   Amperes/sec (2*pi*I0/period), for electric field


a    = 0.00012;       %   conductor diameter in m
M    = 3;          %   number of cross-section subdivisions 
flag = 1;           %   circular cross-section    
sk   = 0;           %   surface current distribution (skin layer)
twist = 0.0;        %   twist, radian per step

[strcoil, check]        = meshwire(Pcenter, a, a, M, flag, sk, twist);  %    Wire model    
[P, t]                  = meshsurface(Pcenter, a, a, M, flag, twist);   %    CAD mesh (optional, slow)    


figure;
bemf1_graphics_coil_CAD(P, t, 1);
view(-140, 40);
%% 

strcoil.I0      = I0;
strcoil.dIdt    = dIdt;
strcoil.P       = P;
strcoil.t       = t;
save('coil', 'strcoil');


coil1 = 'coil_manon.xlsx';
writematrix(strcoil.Pwire,coil1,'Sheet',1)
%% 

mu0                         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)
prec                        = 1e-4;
Inductance                  = bemf6_inductance_neumann_integral(strcoil, mu0, prec)
if check > 3
    disp('Decrease the ratio AvgSegmentLength/AvgSegmentSpacing for the wire mesh. Inductance accuracy cannot be guaranteed.')
end

