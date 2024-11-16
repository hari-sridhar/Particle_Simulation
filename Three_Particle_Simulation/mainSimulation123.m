% Read the Excel file into MATLAB
cases = readtable('Book1.xlsx');
data = readtable('Case_Parameters.xlsx');
% Display the imported data
disp(cases);

% Prompt user to enter a case number
user_case = input('Enter the case number: ');
% Retrieve the corresponding row for the case number
case_data = data(data.CaseNo == user_case, :);

% Step 5: Display or process the parameters for the selected case
if isempty(case_data) 
    disp('Case number not found.');
else
    disp('Parameters for the selected case:');
    disp(case_data);
    R1 = case_data.R1;
    R2 = case_data.R2;
    R3 = case_data.R3;
    X1 = case_data.X1;
    Y1 = case_data.Y1;
    X2 = case_data.X2;
    Y2 = case_data.Y2;
    X3 = case_data.X3;
    Y3 = case_data.Y3;

    simulateParticles(R1, R2, R3, X1, Y1, X2, Y2, X3, Y3);  
    % Further processing can be added here if needed
end



