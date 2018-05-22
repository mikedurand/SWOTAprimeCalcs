function SpecialCases = CheckSpecialCases(Wobs,stdW)

SNR_W_thresh=2;

SNR_W=(var(Wobs)-stdW^2)/stdW^2;

if SNR_W<SNR_W_thresh
    SpecialCases.LowWidthVar=true;
else
    SpecialCases.LowWidthVar=false;
end

return