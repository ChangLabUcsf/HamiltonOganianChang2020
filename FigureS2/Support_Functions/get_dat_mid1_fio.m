function fio = get_dat_mid1_fio(files, coeff1, coeff2, varargin)
%function [fio] = get_dat_mid1_fio(files, coeff1, coeff2, Nbins, Nparts, Nbins_short)
%
%
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
% larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
%
%
% files : cell array of four file names. Each element is a string, and is
% the full path name to a .dat file holding the nonlinearity data.
%
%    Examples:
%
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_1.dat          
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_2.dat          
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_3.dat          
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_4.dat        
%
% coeff1, coeff2 : coefficient vector which tells the polarity of the 
% mid1 and mid2 filters
%
% Npart = 4 : num of test reps in mid code. Should be 4.
%
% Nbins = 15 : num bins used to obtain projection values in mid code. The
% nonlinearity will then have length = 13.
%
% Nbins_short = 14 : num bins used to get nonlinearity values
%
%
% If you call get_data_sta_fio.m as
%
% [fio] = get_dat_mid1_fio(files, coeff1, coeff2);
%
% then Nparts = 4, Nbins = 21, and Nbins_short = 14 by default.
%
% fio is a struct that contains the nonlinearity results. fio has the
% following structure:
%
% fio.nbins = Nbins;
% fio.nbins_short = Nbins_short;
% fio.px1x2_mtx = px1x2_mtx;                 P(x1,x2)
% fio.px1x2spk_mtx = px1x2spk_mtx;           P(x1,x2|spk)
% fio.px1x2spk_pspk_mtx = px1x2spk_pspk_mtx; P(x1,x2|spk) * P(spk)
% fio.pspk_mtx = pspk_mtx;                   P(spk)
% fio.pspkx1x2_mtx = pspkx1x2_mtx;           P(spk|x1,x2)
% fio.x1_mtx = x1_mtx;                       x1
% fio.px1_mtx = px1_mtx;                     P(x1)
% fio.pspkx1_mtx = pspkx1_mtx;               P(spk|x1)
% fio.px1spk_mtx = px1spk_mtx;               P(x1|spk)
% fio.px1spk_pspk_mtx = px1spk_pspk_mtx;     P(x1|spk) * P(spk)
% fio.x1_mean = x1_mean;                     x1
% fio.pspk_mean = pspk_mean;                 P(spk)
% fio.px1_mean = px1_mean;                   P(x1)
% fio.px1_std = px1_std;                     P(x1)
% fio.px1spk_mean = px1spk_mean;             P(x1|spk)
% fio.px1spk_std = px1spk_std;               P(x1|spk)
% fio.pspkx1_mean = pspkx1_mean;             P(spk|x1)
% fio.pspkx1_std = pspkx1_std;               P(spk|x1)
% fio.info;
%
% caa 2/17/09


x1_mtx = [];
pspkx1_mtx = [];
px1_mtx = [];
px1spk_mtx = [];
px1spk_pspk_mtx = [];

pspkx1x2_mtx=[];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];


minpx = 0; % probability can't be lower than 0

if ( isempty(varargin) )
    Nparts = 4;
    Nbins = 15;
    Nbins_short = 14;
elseif ( length(varargin)==1 )
    Nparts = 4;
    Nbins_short = 14;
    Nbins = varargin{1};
elseif ( length(varargin)==2 )
    Nbins = varargin{1};
    Nparts = varargin{2};
    Nbins_short = 14;
elseif ( length(varargin)==3 )
    Nbins = varargin{1};
    Nparts = varargin{2};
    Nbins_short = varargin{3};
elseif ( length(varargin)==4 )
    Nbins = varargin{1};
    Nparts = varargin{2};
    Nbins_short = varargin{3};
    minpx = varargin{4};
end


if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   display('wrong length of coeff vector');
   return
end

if ( length(files) ~= Nparts )
   length(files)
   Nparts
   display('Nparts and num of files do not match');
   return
end


for i = 1:Nparts
   %fname = sprintf('%s_%u.dat',fname_first,i);
   fp = fopen(files{i},'r');
   
   if ( fp == -1 )
      display('error opening file');
      display(files{i});
      return
   end

   x1 = fread(fp, Nbins, 'double'); % first filter projection values
   x2 = fread(fp, Nbins, 'double'); % second filter projection values
   px1x2 = fread(fp, Nbins*Nbins, 'double'); % prior distr
   px1x2spk_pspk = fread(fp, Nbins*Nbins, 'double'); % really rbar * pxt
   pspk = fread(fp, 1, 'double');
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike
   Neff = fread(fp, 1, 'double');

   fclose(fp);

   % Get rid of machine error values
   ind0 = find(px1x2+eps<minpx);
   px1x2spk(ind0) = 0;

   % Reshape into a matrix so that marginals can be calculated
   px1x2 = reshape(px1x2, Nbins, Nbins);  %p(x1,x2)
   px1x2spk = reshape(px1x2spk, Nbins, Nbins);
   px1x2spk_pspk = reshape(px1x2spk_pspk, Nbins, Nbins);


   if ( coeff1(i) == 1 )
       px1x2_r = px1x2;
       px1x2spk_r = px1x2spk;
       x1_r = x1;
   else
       x1_r = -x1;
   end

   if ( coeff2(i) == -1 )
       x2_r = -x2;
   end

   if ( coeff2(i) == 1 )
      x2_r = x2;
   end


   % assign the complete probability distributions
   px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk,Nbins*Nbins,1)]; % P(x1,x2|spike)
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
   pspk_mtx = [pspk_mtx pspk]; % P(spike)
   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)

   % assign the x1 marginal probability distributions

   % P(x1|spike) * P(spike)
   px1spk_pspk = sum(px1x2spk_pspk,2); % sum across columns to get the first marginal
   px1spk_pspk = px1spk_pspk(:);
   px1spk_pspk_mtx = [px1spk_pspk_mtx px1spk_pspk];

   % P(x1|spike)
   px1spk = sum(px1x2spk,2);
   px1spk = px1spk(:);
   px1spk_mtx = [px1spk_mtx px1spk]; 

   % P(x1)
   px1 = sum(px1x2,2);
   px1 = px1(:);
   px1_mtx = [px1_mtx px1]; 

% %    ior1 = px1sp_psp ./ (px1+eps);
% %    ior1 = ior1(:);
% %    ior1_mtx = [ior1_mtx ior1];


   % P(spike|x1)
   pspkx1 = px1spk_pspk ./ (px1+eps);
   pspkx1 = pspkx1(:);
   pspkx1_mtx = [pspkx1_mtx pspkx1]; 

   xm1 = sum(x1_r .* px1); % expected value: x * p(x) of the projection
   x1_r = x1_r - xm1; % center projections on 0
   x1_r = x1_r(:);
   x1_mtx = [x1_mtx x1_r];

end % (for i)


if ( isempty(pspkx1x2_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end


xmin = min(min(x1_mtx));
xmax = max(max(x1_mtx));
maxmax = max([abs(xmin) xmax]);
edges = linspace(-maxmax, maxmax, Nbins_short);

% ior1_rescaled = zeros(Nbins_short-1,Nparts);
pspkx1_rescaled = zeros(Nbins_short-1,Nparts);
px1_rescaled = zeros(Nbins_short-1,Nparts);
px1spk_rescaled = zeros(Nbins_short-1,Nparts);
npoints = zeros(Nbins_short-1,Nparts);

for i = 1:Nbins_short-1
   for j = 1:Nparts
      ind = find( (x1_mtx(:,j)>=edges(i)) & (x1_mtx(:,j)<edges(i+1)) );
      if ( ~isempty(ind) )
         npoints(i,j)=1; %npoints(i,j)+length(ior1_mtx(ind,j));
%          ior1_rescaled(i,j)= ior1_rescaled(i,j) + sum(ior1_mtx(ind,j) .* px1_mtx(ind,j));
         pspkx1_rescaled(i,j)= pspkx1_rescaled(i,j) + sum(pspkx1_mtx(ind,j) .* px1_mtx(ind,j));
         px1_rescaled(i,j)= px1_rescaled(i,j) + sum(px1_mtx(ind,j));
         px1spk_rescaled(i,j) = px1spk_rescaled(i,j) + sum(px1spk_mtx(ind,j));
      end
   end
   x1_mean(i) = 0.5*(edges(i)+edges(i+1));%mean(x1_mtx');
end

nsamples = sum(npoints');

% P(spike|x1)
% ---------------------------------------------------
pspkx1_rescaled = pspkx1_rescaled ./ (px1_rescaled+eps);
pspkx1_mean = sum(pspkx1_rescaled') ./ (nsamples+eps);
pspkx1_std = sqrt( (sum(pspkx1_rescaled'.^2) ./ (nsamples+eps)-pspkx1_mean.^2) ./ (nsamples-1) .* (nsamples+eps) );

% ior1_rescaled = ior1_rescaled./(px1_rescaled+eps);
% ior1_mean = sum(ior1_rescaled')./(nsamples+eps);
% ior1_std = sqrt((sum(ior1_rescaled'.^2)./(nsamples+eps)-ior1_mean.^2)./(nsamples-1).*(nsamples+eps));

% P(x1)
% ---------------------------------------------------
px1_mean = sum(px1_rescaled') ./ (nsamples+eps);
px1_std = sqrt( (sum(px1_rescaled'.^2) ./ (nsamples+eps)-px1_mean.^2) ./ (nsamples-1) );
sumpx1 = sum(px1_mean);
px1_mean = px1_mean ./ sum(px1_mean);

% dx = (max(x1_mean)-min(x1_mean)) / (Nbins_short-1);
% px1_mean = px1_mean ./ dx;
% px1_std = px1_std ./ dx;



% P(x1|sp)
% ---------------------------------------------------
px1spk_mean = sum(px1spk_rescaled') ./ (nsamples+eps);
px1spk_std = sqrt( (sum(px1spk_rescaled'.^2) ./ (nsamples+eps)-px1_mean.^2) ./ (nsamples-1) );
sumpx1spk = sum(px1spk_mean);
px1spk_mean = px1spk_mean ./ sum(px1spk_mean);

pspk_mean = mean(pspk_mtx);



% % Calculate information using Tanya's early method
% % ---------------------------------------------------
% 
% fio = pspkx1_mean / pspk_mean; % F = P(sp|x1) / P(sp) = P(x1|xp) / P(x1)
% 
% index = find(fio>0);
% information = sum( px1_mean(index) .* fio(index) .* log2(fio(index)) )
% 
% index = find(px1spk_mean>0 & px1_mean>0);
% information = sum( px1spk_mean(index) .* log2( px1spk_mean(index) ./ px1_mean(index) ) )
% 
% disp('pause')
% pause




% Calculate information values using a different technique:
% ---------------------------------------------------
pxs1 = px1spk_mtx(:,1);
pxs2 = px1spk_mtx(:,2);
pxs3 = px1spk_mtx(:,3);
pxs4 = px1spk_mtx(:,4);

px1 = px1_mtx(:,1);
px2 = px1_mtx(:,2);
px3 = px1_mtx(:,3);
px4 = px1_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

info = [info1 info2 info3 info4];


% Save output data
% ---------------------------------------------------
fio.nbins = Nbins;
fio.nbins_short = Nbins_short;
fio.px1x2_mtx = px1x2_mtx;
fio.px1x2spk_mtx = px1x2spk_mtx;
fio.px1x2spk_pspk_mtx = px1x2spk_pspk_mtx;
fio.pspk_mtx = pspk_mtx;
fio.pspkx1x2_mtx = pspkx1x2_mtx;
fio.x1_mtx = x1_mtx;
fio.px1_mtx = px1_mtx;
fio.pspkx1_mtx = pspkx1_mtx;
fio.px1spk_mtx = px1spk_mtx;
fio.px1spk_pspk_mtx = px1spk_pspk_mtx;
fio.pspk_mean = pspk_mean;
fio.x1_mean = x1_mean;
fio.px1_mean = px1_mean;
fio.px1_std = px1_std;
fio.px1spk_mean = px1spk_mean;
fio.px1spk_std = px1spk_std;
fio.pspkx1_mean = pspkx1_mean;
fio.pspkx1_std = pspkx1_std;
fio.info = info;

% 
% figure(99); 
% subplot(5,1,3);cla;
% plot(x1,pxs1); hold all;
% plot(x1,px1);
% set(gca,'xlim',[-4, 4]);
% legend('pxs1','px1');
% title(sprintf('MID1 %.3f',mean(info)));

return;



