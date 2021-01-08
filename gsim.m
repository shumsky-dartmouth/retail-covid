
function [Od1i, Od2i, Owi, Sa] = gsim_super(N, L, pi_i, c, ch, cl, r)
% simulation COVID-19 transmission in a grocery store

% N = number of arrivals in the simulation run
% L = length of path
% pi_i = infection rate in population
% c = overeall averade transmission probability, given
% infectious<->suseceptible interactino
% ch = transmission prob of super-spreader
% cl = transmission prob of lower-type
% r = ratio of wake/direct transmission probability (see paper

% Od1i = vector of infections, associated with individual customers, due to
% direct one-way interactions
% Od2i = infectinos associated with two-way direct
% Owi = wake infections
% Sa = arrival times

% Other parameters
pl = 6;
ph = 18;
alpha = 0.7; % probability arrives on the right
pi_n = 0.03;
k1 = 1.15;
k0mult = 2*k1*((pl+ph)/2)*r; % everything needed for k0 except "c"

scale = 1; % number of grocery store 'areas' simulated
           % If scaled to 1, take
           % care of multiple areas by calling multiple times.

% grocery time-period information
dayh = 18; % number of hours in a day
maxlam = 2.23; % peak period arrival rate

% output storage
Od1i = zeros(1,N); % # realized infections due to one-way driect

% opposite direction output storage
Od2i = zeros(1,N); % infections because of opposite direct contacts

Owp = zeros(1,N); % cumulative probability that an arrival is infection via wake
Owi = zeros(1,N); % realized wake infections

% create arrivals using small time-increments; this will allow for changes
% in arrival rate.

% first get the fraction size so that there is a pdel% chance of an arrival in
% the highest-probability time-period
pdel = 0.01;
% now make sure that we have enough time-period to generate N
% multiply by >2 because pdel refers to the max, and we have a triangular
% arrival pattern
Tscale = 2.1;
Tnum = Tscale*N/pdel;
tdel = pdel/maxlam;
Ttimes = (1:Tnum)*tdel; % actual times 

% create vector of arrival probabilities
% setup day-times
daym = dayh*60;
peakt = (2/3)*daym;
slope_up = maxlam/peakt;
slope_down = maxlam/(daym - peakt);

Tprobs = zeros(1,Tnum);
Dtimes = mod(Ttimes,daym);

% generate arrival probs
Tinds = 1:Tnum;
Upinds = Tinds(Dtimes <= peakt);
Tprobs(Upinds) = tdel*(Dtimes(Upinds)*slope_up);
Downinds = Tinds(Dtimes > peakt);
Tprobs(Downinds) = tdel*(maxlam - slope_down*(Dtimes(Downinds)-peakt));

Narrs = (rand(1,Tnum) < Tprobs); % '1' for each arrival
Sa = Ttimes(Narrs); % get arrival times
if (length(Sa) < N)
    error('Not a sufficient number of arrivals generated.');
end;

Sa = Sa(1:N); % take first N of them
Sp = unifrnd(pl,ph,1,N);
Sd = Sa + L ./ Sp;
Sr = rand(1,N) > alpha;
Si = rand(1,N) < pi_i;
Sn = ( (Si==0) & (rand(1,N) < pi_n) );

% now have two levels of 'c' for superspreaders and others.
c_super = ch;
p_super = (c - cl)/(ch - cl);
Sc = ((Si==1) & (rand(1,N) < p_super)) * c_super;
Sc = Sc + ((Sc==0) & (Si==1)) * cl;
% check:
% mean( Sc(Sc>0) );
% length(Sc(Sc>cl))/length(Sc)
% length(Sc(Sc==cl))/length(Sc(Sc>0))
% length(Sc(Sc==c_super))/length(Sc(Sc>0))

% set up for finding interactions
allinds = 1:N;
for (n=1:N)
    if (mod(n,100000)==0)
        sprintf("%d",n)
    end;
    arr = Sa(n);
    dep = Sd(n);
    dir = Sr(n);
    sin = Si(n);
    scn = Sc(n);
    k0n = k0mult*scn;

    % Find the number of direct interactions and infections, same
    % direction. Note that you have to look for infection in either
    % direction, because we are only examining overtakings here (not
    % overtaken) - so each interaction only examined once.
    indsb = ( (Sr==dir) & (Sa < arr) & (Sd > arr) & (Sd > dep) );
    inds = allinds(indsb);
    ninds = length(inds);
    for i=1:ninds
        ind = inds(i);
   %      Od1g(n) = Od1g(n)+1;
   %      Od1n(ind) = Od1n(ind)+1;
        if ( (sin==1) & (Sn(ind)==0) ) % infectious arrival, susceptible customer
            Od1i(ind) = Od1i(ind)+ (rand < scn);
        elseif ( (Sn(n)==0) & (Si(ind)==1) ) % susceptible arrival, infectious customer)
            Od1i(n) = Od1i(n) + (rand < Sc(ind) ) ;
        end;
    end; % end same direction
    
    % Find the number of direct interactions and infections, opposite
    % direction. First, customers already there with an arrival.
    % Note that only have to check if the arrival is infectious, other
    % customer susceptible. The opposite (infectious other customer, susceptible
    % arrival) is automatically checked when the other customer is the
    % arrival; e.g., each interation examined twice.
    indsb = ( (Sr~=dir) & (Sa < arr) & (Sd > arr) );
    inds = allinds(indsb);
    ninds = length(inds);
    for i=1:ninds
        ind = inds(i);
%         Od2a(n) = Od2a(n)+1;
        if ( (sin==1) & (Sn(ind)==0) ) % infectious arrival, susceptible customer
            Od2i(ind) = Od2i(ind)+ (rand < scn);
        end;
    end; 
    
    % Find the number of direct interactions and infections, opposite
    % direction. Second, customers arrive after arrival
    indsb = ( (Sr~=dir) & (Sa > arr) & (Sa < dep) );
    inds = allinds(indsb);
    ninds = length(inds);
    for i=1:ninds
        ind = inds(i);
 %        Od2b(n) = Od2b(n)+1;
        if ( (sin==1) & (Sn(ind)==0) ) % infectious arrival, susceptible customer
            Od2i(ind) = Od2i(ind)+ (rand < scn);
        end;
    end; % end opposite direction
    
    % now wake exposure; can ignore direction, pretend all is one-way
    % like two-way, only have to look at interaction in one direction, because as we
    % run through all arrivals, all interactions are examined twice.
    tq = arr; % to match notation in the notes
    q = Sp(n);
    % category I;
    indsb = ( (Sa < arr) & (Sd > arr) & (Sd < dep) );
    inds = allinds(indsb);
    ninds = length(inds);
    for i=1:ninds
        ind = inds(i);
        if ( (Sn(n)==0) & (Si(ind)==1) ) % susceptible arrival, infectious customer)
            tp = Sa(ind);
            p = Sp(ind);
            k0ind = k0mult*Sc(ind);
            t0 = L/p - (tq - tp);
            const1 = exp(-k1*p*(tq-tp));
            int1 = (k0ind/((p-q)*k1))*(1-exp(-(p-q)*k1*t0));
            Owp1 = const1*int1;
            Owp(n) = Owp(n)+ Owp1;
        end;
    end; % end wake part I
    
    % wake category II;
    indsb = ( (Sa < arr) & (Sd > dep) );
    inds = allinds(indsb);
    ninds = length(inds);
    for i=1:ninds
        ind = inds(i);
        if ( (Sn(n)==0) & (Si(ind)==1) ) % susceptible arrival, infectious customer)
            tp = Sa(ind);
            p = Sp(ind);
            k0ind = k0mult*Sc(ind);
            t0 = p*(tq-tp)/(q-p);
            const1 = exp(-k1*p*(tq-tp));
            int1 = (k0ind/((p-q)*k1))*(1-exp(-(p-q)*k1*t0));
            Owp1 = const1*int1;
            Owp(n) = Owp(n)+ Owp1;
        end;
    end; % end wake part II
    
    % wake category III;
    indsb = ( (Sa > arr) & (Sd < dep) );
    inds = allinds(indsb);
    ninds = length(inds);
    for i=1:ninds
        ind = inds(i);
        if ( (Sn(n)==0) & (Si(ind)==1) ) % susceptible arrival, infectious customer)
            tp = Sa(ind);
            p = Sp(ind);
            k0ind = k0mult*Sc(ind);
            t0 = L/p;
            int1 = (k0ind/((p-q)*k1))*(1-exp(-(p-q)*k1*t0));
            Owp1 = int1;
            Owp(n) = Owp(n)+ Owp1;
        end;
    end; % end wake part III
    
    % generate actual wake infections
    Owi = rand(1,N) < Owp;
end;

end % function

        
    