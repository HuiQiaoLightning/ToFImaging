function bp = Find_Main_Peak_without_Thresh(ap)
[maxvalue,pos] = findpeaks(ap);
[maxvalue,zeros_pos] = max(maxvalue);
pos = pos(zeros_pos);
bp = zeros(size(ap)); 
for n = pos:-1:1
    if ap(n) <= 0
        break;
    else
        bp(n) = ap(n);
    end
end
for n = pos+1:length(ap)
    if ap(n) <= 0
        break;
    else
        bp(n) = ap(n);
    end
end
return