function [PerformanceMeasures] = EvaluatePerformance(BW, GT)

GT = GT(:,:,1);
BW = BW(:,:,1);

if (max(max(GT)) > 40)
    GT = GT > 40;
end

AandB = GT.*BW;
AorB = double((GT+BW) > 0);

notA_and_notB = (1-GT).*(1-BW);
notAandB = (1-GT).*BW;
Aand_notB = GT.*(1-BW);

n_AandB = numel(nonzeros(AandB));
n_notA_and_notB = numel(nonzeros(notA_and_notB));
n_AorB = numel(nonzeros(AorB));
n_notAandB = numel(nonzeros(notAandB));
n_Aand_notB = numel(nonzeros(Aand_notB));

S = n_AandB / n_AorB;
FP = n_notAandB / n_AorB;
FN = n_Aand_notB / n_AorB;

Precision = n_AandB / (n_AandB + n_notAandB + eps);
Recall = n_AandB / (n_AandB + n_Aand_notB + eps);
Accuracy = n_AandB / (n_AorB + eps);
Fmeasure = 2*((Precision*Recall)/(Precision+Recall+eps));
PerformanceMeasures = [S FP FN Precision Recall Accuracy Fmeasure];
end