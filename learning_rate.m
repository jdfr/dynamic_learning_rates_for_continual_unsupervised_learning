function LearningRate = learning_rate(FunctionType, FunctionConstant, consensus_error, NdxStep, NumSteps, MaxLearningRate, MinLearningRate)
    switch FunctionType
      case 1
        LearningRate = (consensus_error)*FunctionConstant;
      case 2
        LearningRate = (consensus_error.^2)*FunctionConstant;
      case 3
        LearningRate = (consensus_error.^3)*FunctionConstant;
      case 4
        LearningRate = exp(consensus_error)*FunctionConstant;
      case 5
        if NdxStep<0.5*NumSteps   
            % Ordering phase: linear decay
            LearningRate=0.4*(1-NdxStep/NumSteps);
        else
            % Convergence phase: constant
            LearningRate=0.01;        
        end
    end
    LearningRate = max(MinLearningRate, min(MaxLearningRate, LearningRate));

