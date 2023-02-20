%% Main code for the decoupling of HMQC by LRD
% Author : Tianyu Qiu
% Email  : tianyuqiu@stu.xmu.edu.cn
% Feb. 16, 2023.
%% 
% This is the matlab scripts folder used in the following paper:
% Tianyu Qiu, Amir Jahangiri, Xiao Han, Dmitry Lesovoy, Tatiana Agback, Peter Agback, Adnane Achour, Xiaobo Qu*, Vladislav Orekhov*
% Resolution enhancement of NMR by decoupling with low-rank Hankel model.
% Corresponding Author: Xiaobo Qu and Vladislav Orekhov
% Email: quxiaobo@xmu.edu.cn, vladislav.orekhov@nmr.gu.se
%% 
% We make our software routines available for non-profit scientific research, enabling others researchers to understand, reproduce and extend our work. 
% All rights are reserved by the authors. Unauthorized use of the routines for industrial or profit-oriented activities is expressively prohibited.
% If you use this Comprehensive codes, please cite the papers below:

%[1] Xiaobo Qu*, Maxim Mayzel, Jian-Feng Cai, Zhong Chen, Vladislav Orekhov*. 
% Accelerated NMR spectroscopy with low-rank reconstruction, Angewandte Chemie International Edition, 54(3):852-854, 2015.

%[2] Tianyu Qiu, Amir Jahangiri, Xiao Han, Dmitry Lesovoy, Tatiana Agback, Peter Agback, Adnane Achour, Xiaobo Qu*, Vladislav Orekhov*. 
% Resolution enhancement of NMR by decoupling with low-rank Hankel model, arXiv preprint arXiv:2212.01144, 2022.
%%
clear;
close all;
%%
currentFolder = pwd;
addpath(genpath(currentFolder));
%% import real NMR data
load 2D_JCoup.mat
[row,col] = size(signal);
FID_J = signal./max(abs(signal(:))); % normalization
%% construct vector c
J = 35;   % Assume that J=35 Hz for one-bond couplings occur between adjacent 13C atoms, e.g., C_alpha-C_beta
dt=1/4299.226;  % Read from raw data
t0=0;
tEnd=(col-1)*dt+t0;
t=t0:dt:tEnd;
c = cos(pi*t*J).';
%% set mask
mask = ones(size(col,1)); % for full sampling, mask is a indentity matrix. For convenience, we use a vector whose elements are 1 to replace it.
%% set regularization parameter
lamda = 4e2;
%% decoupled by LRD
for ind = 1:row
    [FID_decoup_slice] = solver_LRD(FID_J(ind,:).',c,mask,FID_J(ind,:).',lamda); 
    FID_decoup(:,ind) = FID_decoup_slice;    
end
%% obtain Fourier spectra & normalization
Spec = flipud((fft(FID_J,[],2)).');
spec_norm = abs(Spec)/max(abs(Spec(:)));
DecoupSpec = flipud(fft(FID_decoup,[],1));
Decoup_norm = abs(DecoupSpec)/max(abs(DecoupSpec(:)));
%% present 2D J-coupled & decoupled spectra
close all
figure(1);
contour(spec_norm,20);title('J-Coupled spectrum');grid on;
figure(2);
contour(Decoup_norm,20);title('Decoupled spectrum');grid on;