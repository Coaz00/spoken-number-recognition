function LPC = feature_extraction(x,fs)
    
    win = 20*1e-3*fs;
    p = 15;

    num = round(length(x)/win);
    LPC = zeros(p + 1, num);

    k = 1;
    for i = 1:win:(length(x)-win)
        [a, ~] = autocorel(x(i:i+win-1),p); % odredjivanje parametara AR modela autokorelacionom metodom
        LPC(:, k) = transpose(a);
        k = k+1;
    end
end