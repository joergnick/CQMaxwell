clear all
load('data/Err_data_delta10.mat')
%load('data/Err_data_delta0p1.mat')
%load('data/Ref_Err_plot_data.mat')
%Plot_error_matrices_final( Err_H1,Err_L2,h_s,dof_s/2,tau_s,"Absorbing b.c.",2,1)

ERR2=ERR;
load('data/Err_data_delta01.mat')
%load('data/Err_data_newsigns2.mat')

MaxwellErrors(ERR,ERR2,h_s,tau_s, "Absorbing b.c.", 3,1)
