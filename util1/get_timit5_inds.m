function [timit5_inds, timit5] = get_timit5_inds(out)
if nargin<1
    out = [];
end

timit5 = {'fcaj0_si1479', 'fcaj0_si1804', 'fdfb0_si1948', 'fdxw0_si2141', 'fisb0_si2209', 'mbbr0_si2315', 'mdlc2_si2244', 'mdls0_si998', 'mjdh0_si1984', 'mjmm0_si625'};

if ~isempty(out)
    [i,timit5_inds]=intersect({out.name}, timit5);
else
    timit5_inds = [];
end