function [freq,wavl] = frequencies

c = 299792458;
freq = zeros(rinex.max_sys_index,2);
wavl = zeros(rinex.max_sys_index,2);
glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0 0 0 0 0 0 0 0 0 0 0];

for i=1:105
    if i<=rinex.max_sv_index('gps')     % GPS
        freq(i,1) = 10.23*10^6*154;    %Hz
        wavl(i,1) = c/(10.23*10^6*154);%m
        freq(i,2) = 10.23*10^6*120;    %Hz
        wavl(i,2) = c/(10.23*10^6*120);%m
    elseif i<=rinex.max_sv_index('glonass') % GLONASS
        freq(i,1) = (1602 + 0.5625*glok(i-rinex.max_sv_index('gps')))*10^6;    %Hz
        wavl(i,1) = c/((1602 + 0.5625*glok(i-rinex.max_sv_index('gps')))*10^6);%m
        freq(i,2) = (1246 + 0.4375*glok(i-rinex.max_sv_index('gps')))*10^6;    %Hz
        wavl(i,2) = c/((1246 + 0.4375*glok(i-rinex.max_sv_index('gps')))*10^6);%m
    elseif i<=rinex.max_sv_index('galileo') % GALILEO
        freq(i,1) = 10.23*10^6*154;    %Hz
        wavl(i,1) = c/(10.23*10^6*154);%m
        freq(i,2) = 10.23*10^6*115;    %Hz
        wavl(i,2) = c/(10.23*10^6*115);%m
    else        % BEIDOU
        freq(i,1) = 10.23*10^6*152.6;    %Hz
        wavl(i,1) = c/(10.23*10^6*152.6);%m
        freq(i,2) = 10.23*10^6*118;
        wavl(i,2) = c/(10.23*10^6*118);%m
    end
end
end

