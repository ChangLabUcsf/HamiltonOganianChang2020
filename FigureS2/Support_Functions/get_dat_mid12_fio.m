function [fio] = get_dat_mid12_fio(files, coeff1, coeff2, varargin)
%
%function [fio] = get_dat_mid12_fio(files, coeff1, coeff2, Nbins, Nparts)
%
%Get 2D MID nonlinearity
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
% Npart = 4 : num of test reps in mid code. Should always be 4.
%
% Nbins = 15 : num bins used to obtain projection values in mid code. For
% the ICC analysis Nbins = 21.
%
% If you call get_data_sta_fio.m as
%
% [fio] = get_dat_mid1_fio(files, coeff1, coeff2);
%
% then Nparts = 4 and Nbins = 15, by default.
%
% fio is a struct that contains the nonlinearity results. fio has the
% following structure:
%
% fio.nbins = Nbins;
% fio.pspk_mtx = pspk_mtx;                   P(spk)
% fio.x1_mtx = x1_mtx;                       x1
% fio.x2_mtx = x2_mtx;                       x2
% fio.px1x2_mtx = px1x2_mtx;                 P(x1,x2)
% fio.px1x2spk_mtx = px1x2spk_mtx;           P(x1,x2|spk)
% fio.px1x2spk_pspk_mtx = px1x2spk_pspk_mtx; P(x1,x2|spk) * P(spk)
% fio.pspkx1x2_mtx = pspkx1x2_mtx;           P(spk|x1,x2)
% fio.pspk_mean = pspk_mean;                 P(spk)
% fio.x1_mean = x1_mean;                     x1
% fio.x2_mean = x1_mean;                     x2
% fio.px1x2_mean = px1x2_mean;               P(x1,x2)
% fio.px1x2_std = px1x2_std;                 P(x1,x2)
% fio.px1x2spk_mean = px1x2spk_mean;         P(x1,x2|spk)
% fio.px1x2spk_std = px1x2spk_std;           P(x1,x2|spk)
% fio.pspkx1x2_mean = pspkx1x2_mean;         P(spk|x1,x2)
% fio.pspkx1x2_std = pspkx1x2_std;           P(spk|x1,x2)
% fio.info;
%
% caa 2/17/09


x1_mtx = [];
x2_mtx = [];

pspkx1x2_mtx = [];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];

if isempty(varargin)
   Nparts = 4;
   Nbins = 15
elseif ( length(varargin) == 1 )
   Nbins = varargin{1};
   Nparts = 4;
elseif ( length(varargin) == 2 )
   Nbins = varargin{1};
   Nparts = varargin{2};
end


if (length(coeff1)~=Nparts)
    length(coeff1)
    display('wrong length of coeff vector');
    return
end


if ( length(files) ~= Nparts )
   length(files)
   Nparts
   display('Nparts and num of files do not match');
   return
end




for i=1:Nparts

    fp = fopen(files{i},'r');

    if (fp==-1)
        error(sprintf('error opening file %s', files{i}));
        return
    end

   x1 = fread(fp,Nbins,'double');
   x2 = fread(fp,Nbins,'double');
   px1x2 = fread(fp,Nbins*Nbins,'double');
   px1x2spk_pspk = fread(fp,Nbins*Nbins,'double');
   pspk = fread(fp, 1, 'double');

   px1x2 = reshape(px1x2,Nbins,Nbins); % make it a matrix
   px1x2spk_pspk = reshape(px1x2spk_pspk,Nbins,Nbins); % make it a matrix
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike

   fclose(fp);

    if (coeff1(i)==1)
        px1x2_r = px1x2;
        px1x2spk_r = px1x2spk;
        px1x2spk_pspk_r = px1x2spk_pspk;
        x1_r = x1;
    else
        for j=1:Nbins
            px1x2_r(Nbins+1-j,:) = px1x2(j,:);
            px1x2spk_pspk_r(Nbins+1-j,:) = px1x2spk_pspk(j,:);
            px1x2spk_r(Nbins+1-j,:) = px1x2spk(j,:);
            x1_r = -x1;
        end
    end

   if (coeff2(i)==-1)
      for j=1:Nbins
         px1x2_r(:,Nbins+1-j) = px1x2(:,j);
         px1x2spk_pspk_r(:,Nbins+1-j) = px1x2spk_pspk(:,j);
         px1x2spk_r(:,Nbins+1-j) = px1x2spk(:,j);
         x2_r = -x2;
      end
   end

   if (coeff2(i)==1) x2_r = x2; end

   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk_r./(px1x2_r+eps),Nbins*Nbins,1)];

   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk_r,Nbins*Nbins,1)]; % make it a column vector
   px1x2_mtx = [px1x2_mtx reshape(px1x2_r,Nbins*Nbins,1)]; % make it a column vector
   pspk_mtx = [pspk_mtx pspk];
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk_r,Nbins*Nbins,1)]; % make it a column vector

   x1_r = reshape(x1_r,Nbins,1); % stolbets
   x2_r = reshape(x2_r,Nbins,1); % stolbets

   x1_mtx = [x1_mtx x1_r];
   x2_mtx = [x2_mtx x2_r];

end


if (isempty(px1x2spk_pspk_mtx) )
   error('empty px1x2spk_pspk_mtx');
end

xmin = min(min(x1_mtx));
xmax = max(max(x1_mtx));
maxmax = max([abs(xmin) xmax]);
xedges = linspace(-maxmax, maxmax, Nbins+1);

ymin = min(min(x2_mtx));
ymax = max(max(x2_mtx));
maxmax = max([abs(ymin) ymax]);
yedges = linspace(-maxmax, maxmax, Nbins+1);

pspkx1x2_rescaled = zeros(Nbins*Nbins,Nparts);
px1x2_rescaled = zeros(Nbins*Nbins,Nparts);
px1x2spk_rescaled = zeros(Nbins*Nbins,Nparts);

for i=1:Nbins%-1
    for k=1:Nbins%-1
        for j=1:Nparts
            xind = find((x1_mtx(:,j)>=xedges(i))&(x1_mtx(:,j)<xedges(i+1)));
            yind = find((x2_mtx(:,j)>=yedges(k))&(x2_mtx(:,j)<yedges(k+1)));
            temp = 0;
            sum_px1x2 = 0;
            sum_px1x2spk = 0;
            for ii=1:length(xind)
                for jj=1:length(yind)
                    temp = temp + pspkx1x2_mtx((yind(jj)-1)*Nbins+xind(ii),j) .* px1x2_mtx((yind(jj)-1)*Nbins+xind(ii),j);
                    sum_px1x2 = sum_px1x2 + sum(px1x2_mtx((yind(jj)-1)*Nbins+xind(ii),j));
                    sum_px1x2spk = sum_px1x2spk + sum(px1x2spk_mtx((yind(jj)-1)*Nbins+xind(ii),j));
                end
            end
            if (temp>0)
                temp = temp/sum_px1x2;
            end
            pspkx1x2_rescaled((k-1)*Nbins+i,j) = pspkx1x2_rescaled((k-1)*Nbins+i,j) + temp;
            px1x2_rescaled((k-1)*Nbins+i,j) = px1x2_rescaled((k-1)*Nbins+i,j) + sum_px1x2;
            px1x2spk_rescaled((k-1)*Nbins+i,j) = px1x2spk_rescaled((k-1)*Nbins+i,j) + sum_px1x2spk;
        end
    end

    x1_mean(i) = 0.5*(xedges(i)+xedges(i+1));
    x2_mean(i) = 0.5*(yedges(i)+yedges(i+1));

end


% for i = 1:Nbins_short-1
%    for j = 1:Nparts
%       ind = find( (x1_mtx(:,j)>=edges(i)) & (x1_mtx(:,j)<edges(i+1)) );
%       if ( ~isempty(ind) )
%          npoints(i,j)=1; %npoints(i,j)+length(ior1_mtx(ind,j));
% %          ior1_rescaled(i,j)= ior1_rescaled(i,j) + sum(ior1_mtx(ind,j) .* px1_mtx(ind,j));
%          pspkx1_rescaled(i,j)= pspkx1_rescaled(i,j) + sum(pspkx1_mtx(ind,j) .* px1_mtx(ind,j));
%          px1_rescaled(i,j)= px1_rescaled(i,j) + sum(px1_mtx(ind,j));
%          px1spk_rescaled(i,j) = px1spk_rescaled(i,j) + sum(px1spk_mtx(ind,j));
%       end
%    end
%    x1_mean(i) = 0.5*(edges(i)+edges(i+1));
% end


% P(spike|x1,x2)
% ---------------------------------------------------
pspkx1x2_mean = mean(pspkx1x2_rescaled');
pspkx1x2_mean = reshape(pspkx1x2_mean, Nbins, Nbins)';
pspkx1x2_mean = pspkx1x2_mean ./ sum(sum(pspkx1x2_mean));

pspkx1x2_std = sqrt(var(pspkx1x2_rescaled')/Nparts);
pspkx1x2_std = reshape(pspkx1x2_std, Nbins, Nbins)';


% P(x1,x2)
% ---------------------------------------------------
px1x2_mean = mean(px1x2_rescaled');
px1x2_mean = reshape(px1x2_mean, Nbins, Nbins)';
px1x2_mean = px1x2_mean ./ sum(sum(px1x2_mean));

px1x2_std = sqrt(var(px1x2_rescaled')/Nparts);
px1x2_std = reshape(px1x2_std, Nbins, Nbins)';


% P(x1,x2|spike)
px1x2spk_mean = mean(px1x2spk_rescaled');
px1x2spk_mean = reshape(px1x2spk_mean, Nbins, Nbins)';
px1x2spk_mean = px1x2spk_mean ./ sum(sum(px1x2spk_mean));

px1x2spk_std = sqrt(var(px1x2spk_rescaled')/Nparts);
px1x2spk_std = reshape(px1x2spk_std, Nbins, Nbins)';

x12 = x1_mean;
y12 = x2_mean;

pspk_mean = mean(pspk_mtx);


% F(x1,x2) = P(spike|x1,x2) / P(spike) = P(x1,x2|spike) / P(x1,x2)
% ---------------------------------------------------
fx1x2 = pspkx1x2_mean ./ pspk_mean; 
fx1x2 = px1x2spk_mean ./ px1x2_mean; 

index = find(fx1x2>0);
information = sum( sum( px1x2_mean(index) .* fx1x2(index) .* log2(fx1x2(index)) ) );


% Calculate information values using a different technique:
% ---------------------------------------------------
pxs1 = px1x2spk_mtx(:,1);
pxs2 = px1x2spk_mtx(:,2);
pxs3 = px1x2spk_mtx(:,3);
pxs4 = px1x2spk_mtx(:,4);

px1 = px1x2_mtx(:,1);
px2 = px1x2_mtx(:,2);
px3 = px1x2_mtx(:,3);
px4 = px1x2_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

info = [info1 info2 info3 info4];


fio.nbins = Nbins;
fio.pspk_mtx = pspk_mtx;
fio.x1_mtx = x1_mtx;
fio.x2_mtx = x2_mtx;
fio.px1x2_mtx = px1x2_mtx;
fio.px1x2spk_mtx = px1x2spk_mtx;
fio.px1x2spk_pspk_mtx = px1x2spk_pspk_mtx;
fio.pspkx1x2_mtx = pspkx1x2_mtx;
fio.pspk_mean = pspk_mean;
fio.x1_mean = x1_mean;
fio.x2_mean = x2_mean;
fio.px1x2_mean = px1x2_mean;
fio.px1x2_std = px1x2_std;
fio.px1x2spk_mean = px1x2spk_mean;
fio.px1x2spk_std = px1x2spk_std;
fio.pspkx1x2_mean = pspkx1x2_mean;
fio.pspkx1x2_std = pspkx1x2_std;
fio.info = info;


% 
% figure(99); 
% subplot(5,1,5);cla;
% imagesc(reshape(pxs1-px1, [Nbins, Nbins]));
% title(sprintf('MID1+2 %.3f',mean(info)));

return;


