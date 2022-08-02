

function [Temp,Var] = CalculateTempAndVariance(T)
%
%
%
%
%
%% Calculator
for CalcType = {'Value','Smoothed'}
    CT = CalcType{1,1};
    for Info = {'Temp','Var'}
        In = Info{1,1};
        % Pre-allocating storage arrays
        A.(In).(CT) = zeros(size(T{1,1}.(CT))); 
        MaskSum = zeros(size(T{1,1}.(CT)));
        StepStore = cell(size(T,1),1);
        % Summing all the bootstrap steps together
        for m=1:1:size(T,1)
            % Calculating the parameter of interest
            Step = StepCalc(T{m,1}.(CT),T{m,2}.(CT),In);
            Mask = ~isnan(Step);
            % Removing nans from the summation
            Step(~Mask) = 0;
            % Summing steps 
            A.(In).(CT) = A.(In).(CT) + Step;
            MaskSum     = MaskSum + Mask;
            % Storing steps
            StepStore{m,1} = A.(In).(CT);
            StepStore{m,2} = MaskSum;
        end
        % Calculating the final product
        A.(In).(CT) = SumCalc(A.(In).(CT),MaskSum,In);
        % If the variable is variance, track its progression
        if strcmp(In,'Var')
            ValOld = nan.*zeros(size(StepStore{1}));
            for m=1:1:size(T,1)
                ValNew = SumCalc(StepStore{m,1},StepStore{m,2},In);
                ValNew(ValNew==inf) = nan;
                % Calculate max change in variance for each iteration
                A.Var.([CT,'MaxChange'])(1,m) = max(max(abs(ValNew - ValOld)));
                % Updating variance
                ValNew(isnan(ValNew)) = 0; ValOld = ValNew;
            end
        end
    end
end 


%% Parsing out data for return
Temp = A.Temp;
Var  = A.Var;
end

function [Val] = StepCalc(A,B,In)
%
% Input: A:  Array containing temperature profile A
%        B:  Array containing temperature profile B
%        In: String (either 'Temp' or 'Var') that defines the calculation
%            of interest
%
% Output: Val: Array of calculated temperature means or variances
%
%% Determining how to calculate the parameter of interest
if strcmp(In,'Temp') 
    Val = (A+B)./2;  % Temperature calulcation (mean of 2 profiles)
else                 
    Val = (A-B).^2;  % Variance calculations
end
end

function [Val] = SumCalc(A,M,In)
%
% Input: A:  Array containing summed temp or variance
%        M:  Array containing mask sums
%        In: String (either 'Temp' or 'Var') that defines the calculation
%            of interest
%
% Output: Val: Array of calculated temperature means or variances
%
%% Determining how to calculate the parameter of interest
if strcmp(In,'Temp') 
    Val = A./M;         % Temperature calulcation (weigmean)
else   
    M(M<2) = 1;         % Making sure anything with 1 value gets removed & not turned negative
    Val = A./2./(M-1);  % Variance calculations
end
end



