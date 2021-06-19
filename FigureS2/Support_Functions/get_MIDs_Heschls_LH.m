function [mtx_sta mtx1 mtx2 stainf mid1inf mid2inf mid12inf strfinf] = get_MIDs_Heschls_LH(date, site, cell,Nparts,rw,col,opt_level1,opt_level2,subj,paper_data_dir,dStims)
%cell = vector of cell numbers.
%set opt_level1 = 1 and opt_level2 = 2;

Nv = rw;
Nh = col;
cx = 1;
cy = 1;

Nbins = 21;
Nbins_medium = 15;
Nbins_short = 14;
%Nparts = 4;


mid_computed = 1;
%--------------------------------------------------------------------
%   Data Location
%--------------------------------------------------------------------
prefix = [pwd '/'];



%--------------------------------------------------------------------
%   STA 
%--------------------------------------------------------------------
file_sta = sprintf('%srpsta_%u_%u_%u_1x%ux%u_1', prefix, date, site, cell, Nv, Nh);

file_sta

[v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_filter(file_sta, Nh, Nv, Nparts);



%--------------------------------------------------------------------
%   STA Nonlinearities
%--------------------------------------------------------------------
file_sta_fio{1} = sprintf('%srpx1pxpxt_sta_%u_%u_%u_1x%ux%u_1_1.dat', prefix, date, site, cell, Nv, Nh);
file_sta_fio{2} = sprintf('%srpx1pxpxt_sta_%u_%u_%u_1x%ux%u_1_2.dat', prefix, date, site, cell, Nv, Nh);
file_sta_fio{3} = sprintf('%srpx1pxpxt_sta_%u_%u_%u_1x%ux%u_1_3.dat', prefix, date, site, cell, Nv, Nh);
file_sta_fio{4} = sprintf('%srpx1pxpxt_sta_%u_%u_%u_1x%ux%u_1_4.dat', prefix, date, site, cell, Nv, Nh);

[fiosta] = get_dat_sta_fio(file_sta_fio, coeff_sta);
stainf = fiosta.info;


%--------------------------------------------------------------------
%   Ridge STRF information 
%--------------------------------------------------------------------
[fiostrf] = get_strfridge_fio(subj, cell, paper_data_dir, Nbins,dStims);
strfinf = fiostrf.info;


%--------------------------------------------------------------------
%   MID1
%--------------------------------------------------------------------

% rpdtest2_prelim_v1_707_1_1x25x20_1_*1.dat
optlevel = opt_level1;

file_v1 = sprintf('%srpdtest%u_v1_%u_%u_%u_1x%ux%u_%u', prefix, optlevel, date, site, cell, Nv, Nh, cy);

[v1, coeff1, projection1, mtx1] = get_auditory_filter(file_v1, Nh, Nv, Nparts);


%--------------------------------------------------------------------
%   MID2
%--------------------------------------------------------------------

% rpdtest2_prelim_v2_707_1_1x25x20_1_*1.dat
optlevel = opt_level2;

file_v2 = sprintf('%srpdtest%u_v2_%u_%u_%u_1x%ux%u_%u', prefix, optlevel, date, site, cell, Nv, Nh, cy);

[v2, coeff2, projection2, mtx2] = get_auditory_filter(file_v2, Nh, Nv, Nparts);



%--------------------------------------------------------------------
%   MID1, MID2, MID12 Nonlinearities
%--------------------------------------------------------------------
files{1} = sprintf('%srpdx1x2px_pxt_1_%u_%u_%u_1x%ux%u_1_1.dat', prefix, date, site, cell, Nv, Nh);
files{2} = sprintf('%srpdx1x2px_pxt_1_%u_%u_%u_1x%ux%u_1_2.dat', prefix, date, site, cell, Nv, Nh);
files{3} = sprintf('%srpdx1x2px_pxt_1_%u_%u_%u_1x%ux%u_1_3.dat', prefix, date, site, cell, Nv, Nh);
files{4} = sprintf('%srpdx1x2px_pxt_1_%u_%u_%u_1x%ux%u_1_4.dat', prefix, date, site, cell, Nv, Nh);

% files{1} = sprintf('%srpdx1x2px_pxt_2_%u_%u_%u_1x%ux%u_1_1.dat', prefix, date, site, cell, Nv, Nh);
% files{2} = sprintf('%srpdx1x2px_pxt_2_%u_%u_%u_1x%ux%u_1_2.dat', prefix, date, site, cell, Nv, Nh);
% files{3} = sprintf('%srpdx1x2px_pxt_2_%u_%u_%u_1x%ux%u_1_3.dat', prefix, date, site, cell, Nv, Nh);
% files{4} = sprintf('%srpdx1x2px_pxt_2_%u_%u_%u_1x%ux%u_1_4.dat', prefix, date, site, cell, Nv, Nh);

Nbins = 21;
Nparts = 4;

[fio1] = get_dat_mid1_fio(files, coeff1, coeff2, Nbins, Nparts);
mid1inf = fio1.info;

[fio2] = get_dat_mid2_fio(files, coeff1, coeff2, Nbins, Nparts);
mid2inf = fio2.info;

[fio12] = get_dat_mid12_fio(files, coeff1, coeff2, Nbins, Nparts);
mid12inf = fio12.info;






