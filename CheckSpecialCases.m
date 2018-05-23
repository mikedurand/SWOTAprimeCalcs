function SpecialCases = CheckSpecialCases(Wobs,stdW,Hobs,stdH)

SNR_W_thresh=1;
SNR_H_thresh=3;

SNR_W=(var(Wobs)-stdW^2)/stdW^2;

if SNR_W<SNR_W_thresh
    SpecialCases.LowWidthVar=true;
else
    SpecialCases.LowWidthVar=false;
end

SNR_H=(var(Hobs)-stdH^2)/stdH^2;

if SNR_H<SNR_H_thresh
    SpecialCases.LowHeightVar=true;
else
    SpecialCases.LowHeightVar=false;
end

if ~SpecialCases.LowWidthVar && ~SpecialCases.LowHeightVar
    SpecialCases.Nominal=true;
else
    SpecialCases.Nominal=false;
end

return