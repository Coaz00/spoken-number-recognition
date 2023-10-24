function [m1, m2, m3, m4, m5, m6] = sekvence(x)

N = length(x);
m1 = zeros(1,N);
m2 = zeros(1,N);
m3 = zeros(1,N);
m4 = zeros(1,N);
m5 = zeros(1,N);
m6 = zeros(1,N);

maxs = zeros(1,N); % sadrzi nule svugde sem na mestima lok. maks.
mins = zeros(1,N); % isto za min.

[peaks, max_idxs] = findpeaks(x);
maxs(max_idxs) = peaks;
max_idxs = max_idxs';

[peaks, min_idxs] = findpeaks(-x); % minimum pronalazimo kao maksimume -x
mins(min_idxs) = -peaks;
min_idxs = min_idxs';

m1 = maxs;
m1(m1<0) = 0;

m4 = -mins;
m4(m4<0) = 0;

maxp = 0; % vrednost prethodnog maksimuma
for i = max_idxs
    if (isempty(find(min_idxs < i, 1, 'last')))
        minp = 0;
    else
        idx = find(min_idxs < i, 1, 'last');
        minp = mins(min_idxs(idx));
    end
    
    m2(i) = max(0, maxs(i)-minp);
    m3(i) = max(0, maxs(i)-maxp);
    maxp = maxs(i);
end

minp = 0; % vrednost prethodnog minimuma
for i = min_idxs
    if (isempty(find(max_idxs < i, 1, 'last')))
        maxp = 0;
    else
        idx = find(max_idxs < i, 1, 'last');
        maxp = mins(max_idxs(idx));
    end
    
    m5(i) = max(0, -(mins(i)-minp));
    m6(i) = max(0, -(mins(i)-maxp));
    minp = mins(i);
end

end