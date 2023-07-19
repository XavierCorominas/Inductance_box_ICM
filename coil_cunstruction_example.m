
% Add Paths
close all 
clear all 
clc

addpath (<<Path to your data coil folder>>) 
addpath (<<Path to your Simnibs4.0 folder>>)



% >reate cells representing coil coordinates and correct for empty values
T = readtable('coilCoordinates1mm.xlsx');
H= table2array(T);

G = H;
idf = isnan(G);
find(idf);
G(idf) = [];
G = reshape(G,[],3);

%Save coil in .xlsx format
cd (<<Path to your data coil folder>>)
coil = 'custom_coil.xlsx';
writematrix(G,coil,'Sheet',1,'Range','D1')



%% Fast visualization of magnetic vectors.
% Once the magetic vercotrs are calculated with the Box ans exported in
% NIFTI format (.nii) we can run the following visualization.

app.hdr = nifti_load('nifti_blank.nii');
            
            app.hdr.sform(3,4)=-85;

            [app.xx,app.yy,app.zz]=meshgrid(linspace(-((app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4)),(app.hdr.dim(2)-1)*app.hdr.sform(1,1)+app.hdr.sform(1,4),app.hdr.dim(2)),...
                                            linspace(-((app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4)),(app.hdr.dim(3)-1)*app.hdr.sform(2,2)+app.hdr.sform(2,4),app.hdr.dim(3)),...
                                            linspace((app.hdr.dim(4)-1)*app.hdr.sform(3,3)+app.hdr.sform(3,4),app.hdr.sform(3,4),app.hdr.dim(4)));


hdr = nifti_load('custom_coil.nii')

                   n=5; %set resolution
                
q=quiver3(app.xx(1:n:end,1:n:end,1:n:end),app.yy(1:n:end,1:n:end,1:n:end),app.zz(1:n:end,1:n:end,1:n:end),hdr.vol(1:n:end,1:n:end,1:n:end,2),hdr.vol(1:n:end,1:n:end,1:n:end,1),hdr.vol(1:n:end,1:n:end,1:n:end,3));
        