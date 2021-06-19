function [v_mean, coeff, projection, mtx] = get_auditory_filter(fname_first, Nh, Nv, varargin)
%
%[v_mean, coeff, projection, mtx] = get_auditory_filter(fname_first, Nh, Nv, varargin)
%
%fname_first : dat filter filename. 
%              Examples: rpsta_707_1_1x25x20_1_1.dat 
%                     or rpdtest2_prelim_v1_707_1_1x25x20_1_4.dat
%
% Nh : number of time bins, or horizontal bins, in filter. This is always
% equal to 20.
%
%Nv : number of vertical bins in filter
%
%varargin : number of test reps. Should be 4.
%
% caa 2/16/09


% These really don't change, ever.
%Nv = 25;
%Nh = 20;
nlags = 1;

if ( ~isempty(varargin) )
    Nparts = varargin{1};
else 
    Nparts = 4;
end

fsize = Nv * Nh;
Nn = fsize * nlags; % total number of elements in the STRF

mtx = [];

for i = 1:Nparts
	fname = sprintf('%s_%u.dat', fname_first, i);
   fp = fopen(fname,'rb');
	if ~(fp== -1)
		[v] = fread(fp,'double');
		v = reshape(v,Nn,1);
		mtx = [mtx,v];
		fclose(fp);
	else
		v = [];
	end
    
end


if (isempty(mtx) )
    fname
    error('empty mtx in plot_a_vector');
end

mtx = reshape(mtx, Nn, Nparts);
coeff(1) = 1;

for i = 2:Nparts
    coeff(i)=sign(sum(mtx(:,i).*mtx(:,1)));
    if ( coeff(i) == -1 )
        mtx(:,i)=mtx(:,i)*coeff(i);
    end
end

v_mean = mean(mtx')';
v_mean = v_mean ./ sqrt(sum(sum(v_mean.*v_mean)));
v_std = sqrt(sum(var(mtx'))*(Nparts-1)/Nn);
v_mean = v_mean ./ v_std;
v_mean = reshape(v_mean,fsize,nlags);

% for testrep=1:Nparts
%     temp_vector=mtx(:,testrep);
%     temp_vector=temp_vector/sqrt(sum(temp_vector.*temp_vector));
%     temp_vector=Nparts*v_mean-(Nparts-1)*temp_vector;%ATTENTION
%     mtx(:,testrep)=temp_vector/sqrt(sum(temp_vector.*temp_vector));
%     if (testrep >1)
%         mtx(:,testrep)=mtx(:,testrep)*sign(sum(mtx(:,testrep).*mtx(:,1)));
%     end
%     
% end
% 
% v_mean=mean(mtx');
% v_mean=v_mean./sqrt(sum(sum(v_mean.*v_mean)));

% if ( nlag_end > nlags )
%     for i=1:Nv*Nh
%         v_mean(i,nlag_end)=0;
%     end    
% end

cm = max([abs(min(min(v_mean))),max(max(v_mean))]);

vprocess = reshape(v_mean(:,1),Nv,Nh);


k = 1;
for i = 1:Nparts
    for j = i+1:Nparts
        projection(k) = sum(mtx(:,i).*mtx(:,j));
        k = k+1;
    end
end

return;
