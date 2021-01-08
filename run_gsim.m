% run simulation functions.  Note that with 400,000 customers each run
% below can require a long time (possibly hours).

% This run has all infectious customers with the same probability "c"
%  of transmission
N=400000;
tic
% set seed
rng(1);
% gsim_super2(N,L, pi_i, c, ch, cl, r)
c_super = 0.001;
sprintf('C-high = %0.4f', c_super)
[O11,O21,Ow1,Sa1] = gsim(N,80,0.006,0.001,c_super,0.0001,0.25);
save samec
[O12,O22,Ow2,Sa2] = gsim(N,80,0.006,0.001,c_super,0.0001,0.25);
save samec
[O13,O23,Ow3,Sa3] = gsim(N,0.7*80,0.006,0.001,c_super,0.0001,0.25);
save samec
if (max([O11,O21,Ow1,O12,O22,Ow2,O13,O23,Ow3])>1)
    sprintf("Repeated infection detected.")
end;
toc
clear

% This run has a 'super-spreader' with a 10% chance of transmission, but
% the same overall value of 'c'
tic
% set seed
rng(1);
% gsim_super2(N,L, pi_i, c, ch, cl, r)
N=400000;
c_super = 0.1;
sprintf('C-high = %0.3f', c_super)
[O11,O21,Ow1,Sa1] = gsim(N,80,0.006,0.001,c_super,0,0.25);
save superspreader
[O12,O22,Ow2,Sa2] = gsim(N,80,0.006,0.001,c_super,0,0.25);
save superspreader
[O13,O23,Ow3,Sa3] = gsim(N,80,0.006,0.001,c_super,0,0.25);
save superspreader
if (max([O11,O21,Ow1,O12,O22,Ow2,O13,O23,Ow3])>1)
    sprintf("Repeated infection detected.")
end;
toc

% look at results; first no super-spreader...
clear;
load samec

% look at summary of results
dayh = 18;
daym = dayh*60;

% combine the infection times in a single stream
mint=min([Sa1(length(Sa1)),Sa2(length(Sa2)),Sa3(length(Sa3))]);
% take only infections within that time period
O11 = logical(O11(Sa1<=mint));
O21 = logical(O21(Sa1<=mint));
O12 = logical(O12(Sa2<=mint));
O22 = logical(O22(Sa2<=mint));
O13 = logical(O13(Sa3<=mint));
O23 = logical(O23(Sa3<=mint));

Ow1 = logical(Ow1(Sa1<=mint));
Ow2 = logical(Ow2(Sa1<=mint));
Ow3 = logical(Ow3(Sa1<=mint));

Sa1t = Sa1(Sa1<=mint);
Sa2t = Sa2(Sa2<=mint);
Sa3t = Sa3(Sa3<=mint);

ndays = mint/daym;
% combine all infections; direct first
SDinf = sort( [Sa1t(O11), Sa1t(O21), Sa2t(O12), Sa2t(O22), Sa3t(O13), Sa3t(O23)] );
sprintf('total direct: %0.2f', length(SDinf)/ndays)
% wake
SWinf = sort( [Sa1t(Ow1), Sa2t(Ow2), Sa3t(Ow3)] );
sprintf('total wake: %0.2f', length(SWinf)/ndays)
sprintf('total: %0.2f', (length(SDinf)+length(SWinf))/ndays)

% copare with base case from Maple, 
% direct = 3.61968/day; wake = 1.8892/day

% get all infections
Sinf = sort( [SDinf,SWinf] );

% create histogram of #/day
ndaysf = floor(mint/daym);
infday = zeros(1,ndaysf);
for (i=1:ndaysf)
    dayl=(i-1)*daym;
    dayh=i*daym;
    infday(i) = sum( (Sinf>=dayl) & (Sinf < dayh) );
end;
maxgrp = 4;
hist = histogram(infday,'Normalization','Probability')
xlabel('Number of infections in a day')
ylabel('Fraction of days')
xlim([-0.75,maxgrp+.75]);
xticks(0:maxgrp);

% plot number in group
nums = hist.Values .* [0:maxgrp];
pnums = nums/sum(nums);
bout = bar(1:maxgrp,[pnums(2:maxgrp+1)],'BarWidth',1);
xlabel('Number of infections in a day (x)')
ylabel('Fraction of infections on a day with x infections')
bout(1).FaceColor = [0.5843,0.8157,0.9882];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% now super-spreader
clear;
load superspreader

% look at summary of results
dayh = 18;
daym = dayh*60;

% combine the infection times in a single stream
mint=min([Sa1(length(Sa1)),Sa2(length(Sa2)),Sa3(length(Sa3))]);
% take only infections within that time period
O11 = logical(O11(Sa1<=mint));
O21 = logical(O21(Sa1<=mint));
O12 = logical(O12(Sa2<=mint));
O22 = logical(O22(Sa2<=mint));
O13 = logical(O13(Sa3<=mint));
O23 = logical(O23(Sa3<=mint));

Ow1 = logical(Ow1(Sa1<=mint));
Ow2 = logical(Ow2(Sa1<=mint));
Ow3 = logical(Ow3(Sa1<=mint));

Sa1t = Sa1(Sa1<=mint);
Sa2t = Sa2(Sa2<=mint);
Sa3t = Sa3(Sa3<=mint);

ndays = mint/daym;
% combine all infections; direct first
SDinf = sort( [Sa1t(O11), Sa1t(O21), Sa2t(O12), Sa2t(O22), Sa3t(O13), Sa3t(O23)] );
sprintf('total direct: %0.2f', length(SDinf)/ndays)
% wake
SWinf = sort( [Sa1t(Ow1), Sa2t(Ow2), Sa3t(Ow3)] );
sprintf('total wake: %0.2f', length(SWinf)/ndays)
sprintf('total: %0.2f', (length(SDinf)+length(SWinf))/ndays)

% copare with base case from Maple, 
% direct = 3.61968/day; wake = 1.8892/day

% get all infections
Sinf = sort( [SDinf,SWinf] );

% create histogram of #/day
ndaysf = floor(mint/daym);
infday = zeros(1,ndaysf);
for (i=1:ndaysf)
    dayl=(i-1)*daym;
    dayh=i*daym;
    infday(i) = sum( (Sinf>=dayl) & (Sinf < dayh) );
end;
maxgrp = 8;
hist = histogram(infday,'Normalization','Probability')
xlabel('Number of infections in a day')
ylabel('Fraction of days')
xlim([-0.75,maxgrp+.75]);
xticks(0:maxgrp);

% plot number in group
nums = hist.Values .* [0:maxgrp];
pnums = nums/sum(nums);
bout = bar(1:maxgrp,[pnums(2:maxgrp+1)],'BarWidth',1);
xlabel('Number of infections in a day (x)')
ylabel('Fraction of infections on a day with x infections')
bout(1).FaceColor = [0.5843,0.8157,0.9882];
