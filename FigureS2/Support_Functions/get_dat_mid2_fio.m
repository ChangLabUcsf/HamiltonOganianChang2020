function [fio] = get_dat_mid2_fio(files, coeff1, coeff2, varargin)
%function [fio] = get_dat_mid2_fio(files, coeff1, coeff2, Nbins, Nparts, Nbins_short)
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
% Nbins = 15 : num bins used to obtain projection values in mid code
%
% Nbins_short = 14 : num bins used to get nonlinearity values. The 
% nonlinearity will then have length = 13.
%
%
% If you call get_data_sta_fio.m as
%
% [fio] = get_dat_mid2_fio(files, coeff1, coeff2);
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
% fio.x2_mtx = x2_mtx;                       x2
% fio.px2_mtx = px2_mtx;                     P(x2)
% fio.pspkx2_mtx = pspkx2_mtx;               P(spk|x2)
% fio.px2spk_mtx = px2spk_mtx;               P(x2|spk)
% fio.px2spk_pspk_mtx = px2spk_pspk_mtx;     P(x2|spk) * P(spk)
% fio.x2_mean = x2_mean;                     x2
% fio.pspk_mean = pspk_mean;                 P(spk)
% fio.px2_mean = px2_mean;                   P(x2)
% fio.px2_std = px2_std;                     P(x2)
% fio.px2spk_mean = px2spk_mean;             P(x2|spk)
% fio.px2spk_std = px2spk_std;               P(x2|spk)
% fio.pspkx2_mean = pspkx2_mean;             P(spk|x2)
% fio.pspkx2_std = pspkx2_std;               P(spk|x2)
% fio.info;
%
% caa 2/17/09


x2_mtx = [];
pspkx2_mtx = [];
px2_mtx = [];
px2spk_mtx = [];
px2spk_pspk_mtx = [];

pspkx1x2_mtx = [];
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


   if ( coeff1(i) == -1 )
       x1_r = -x1;
   else
       x1_r = x1;
   end

   if ( coeff2(i) == -1 )
       x2_r = -x2;
   else
       x2_r = x2;
   end

   if (coeff2(i)==1)
      x2_r = x2;
   end


   % assign the complete probability distributions
   px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk,Nbins*Nbins,1)]; % P(x1,x2|spike)
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
   pspk_mtx = [pspk_mtx pspk]; % P(spike)
   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)


   % assign the x2 marginal probability distributions

   % P(x2|spike) * P(spike)
   px2spk_pspk = sum(px1x2spk_pspk,1); % sum across columns to get the first marginal
   px2spk_pspk = px2spk_pspk(:);
   px2spk_pspk_mtx = [px2spk_pspk_mtx px2spk_pspk];


   % P(x2|spike)
   px2spk = sum(px1x2spk,1);
   px2spk = px2spk(:);
   px2spk_mtx = [px2spk_mtx px2spk];


   % P(x2)
   px2 = sum(px1x2,1);
   px2 = px2(:);
   px2_mtx = [px2_mtx px2];


   % P(spike|x2)
   pspkx2 = px2spk_pspk ./ (px2+eps);
   pspkx2 = pspkx2(:);
   pspkx2_mtx = [pspkx2_mtx pspkx2];


   xm2 = sum( x2_r .* px2 ); % expected value: x * p(x) of the projection
   x2_r = x2_r - xm2;
   x2_r = x2_r(:);
   x2_mtx = [x2_mtx x2_r];


end % (for i)



if ( isempty(pspkx1x2_mtx) )
    error('empty pspkx1x2_mtx in nonlinearity file.');
end


xmin = min(min(x2_mtx));
xmax = max(max(x2_mtx));
maxmax = max([abs(xmin) xmax]);
edges = linspace(-maxmax, maxmax, Nbins_short);


pspkx2_rescaled = zeros(Nbins_short-1,Nparts);
px2_rescaled = zeros(Nbins_short-1,Nparts);
px2spk_rescaled = zeros(Nbins_short-1,Nparts);
npoints = zeros(Nbins_short-1,Nparts);


for i = 1:Nbins_short-1
   for j = 1:Nparts
      ind = find((x2_mtx(:,j)>=edges(i))&(x2_mtx(:,j)<edges(i+1)));
      if ~isempty(ind)
         npoints(i,j) = 1; 
         pspkx2_rescaled(i,j) = pspkx2_rescaled(i,j) + sum( pspkx2_mtx(ind,j) .* px2_mtx(ind,j) );
         px2_rescaled(i,j) = px2_rescaled(i,j) + sum( px2_mtx(ind,j) );
         px2spk_rescaled(i,j) = px2spk_rescaled(i,j) + sum( px2spk_mtx(ind,j) );
      end
   end
   x2_mean(i) = 0.5*(edges(i)+edges(i+1));
end

nsamples = sum(npoints');


% P(spike|x2)
% ---------------------------------------------------
pspkx2_rescaled = pspkx2_rescaled ./ (px2_rescaled+eps);
pspkx2_mean = sum(pspkx2_rescaled') ./ (nsamples+eps);
pspkx2_std = sqrt( (sum(pspkx2_rescaled'.^2) ./ (nsamples+eps)-pspkx2_mean.^2) ./ (nsamples-1) .* (nsamples+eps) );


% P(x2)
% ---------------------------------------------------
px2_mean = sum(px2_rescaled') ./ (nsamples+eps);
px2_std = sqrt( (sum(px2_rescaled'.^2) ./ (nsamples+eps)-px2_mean.^2) ./ (nsamples-1));
sumpx2 = sum(px2_mean);
px2_mean = px2_mean ./ sum(px2_mean);


% P(x2|spike)
% ---------------------------------------------------
px2spk_mean = sum(px2spk_rescaled') ./ (nsamples+eps);
px2spk_std = sqrt( (sum(px2spk_rescaled'.^2) ./ (nsamples+eps)-px2_mean.^2) ./ (nsamples-1) );
sumpx2spk = sum(px2spk_mean);
px2spk_mean = px2spk_mean ./ sum(px2spk_mean);

pspk_mean = mean(pspk_mtx);



% % Calculate information using Tanya's early method
% % ---------------------------------------------------
% fio = pspkx2_mean / pspk_mean; % F = P(spk|x2) / P(spk) = P(x2|xp) / P(x2)
% 
% index = find(fio>0);
% information = sum( px2_mean(index) .* fio(index) .* log2(fio(index)) )
% 
% index = find(px2spk_mean>0 & px2_mean>0);
% information = sum( px2spk_mean(index) .* log2( px2spk_mean(index) ./ px2_mean(index) ) )


% Calculate information values using a different technique:
% ---------------------------------------------------
pxs1 = px2spk_mtx(:,1);
pxs2 = px2spk_mtx(:,2);
pxs3 = px2spk_mtx(:,3);
pxs4 = px2spk_mtx(:,4);

px1 = px2_mtx(:,1);
px2 = px2_mtx(:,2);
px3 = px2_mtx(:,3);
px4 = px2_mtx(:,4);

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
fio.x2_mtx = x2_mtx;
fio.px2_mtx = px2_mtx;
fio.pspkx2_mtx = pspkx2_mtx;
fio.px2spk_mtx = px2spk_mtx;
fio.px2spk_pspk_mtx = px2spk_pspk_mtx;
fio.pspk_mean = pspk_mean;
fio.x2_mean = x2_mean;
fio.px2_mean = px2_mean;
fio.px2_std = px2_std;
fio.px2spk_mean = px2spk_mean;
fio.px2spk_std = px2spk_std;
fio.pspkx2_mean = pspkx2_mean;
fio.pspkx2_std = pspkx2_std;
fio.info = info;
% 
% figure(99); 
% subplot(5,1,4);cla;
% plot(x1,pxs1); hold all;
% plot(x1,px1);
% set(gca,'xlim',[-4, 4]);
% legend('pxs1','px1');
% title(sprintf('MID2 %.3f',mean(info)));

return;







