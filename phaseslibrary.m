function phasesproperties=phaseslibrary()

disp('definition of phase library')

% Define solvent names
PhasesNames = {'protein', 'silica', 'calcite', 'Fe3O4', 'CdSe', 'Au'};

% Define properties as column vectors
PhasesNamesRow = {'protein'; 'silica'; 'calcite'; 'Fe3O4'; 'CdSe'; 'Au'};
Density = [1220,	2285,	2711,	5170,	5820,	19320]'; % kg/m³

% Define row names correctly
RowNames = {
    'Phase'
    'Density [kg/m³]', 
    };

% Create the table
phasesproperties = table(PhasesNamesRow,Density);
phasesproperties.Properties.RowNames=PhasesNames;
phasesproperties.Properties.VariableNames=RowNames;
