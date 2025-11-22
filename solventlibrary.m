function solventproperties=solventlibrary()

disp('definition of solvent library')

% Define solvent names
SolventNames = {'water', 'ethanol', 'toluene', 'tce', 'hexadecane'};

% Define properties as column vectors
SolventNamesRow = {'water'; 'ethanol'; 'toluene'; 'tce'; 'hexadecane'};
MW = [18.01528; 46.069; 92.141; 162.82; 226.448]; % g/mol
solventmass = [2.99158E-26; 7.65012E-26; 1.53007E-25; 2.70375E-25; 3.76035E-25]; % Kg
solventradius = [1.65958E-10; 2.45599E-10; 3.0041E-10; 2.94218E-10; 4.20994E-10]; % m
Density = [1000; 789; 862.3; 1622; 770]; % kg/m³
Viscosity = [8.90E-04; 1.07E-03; 5.60E-04; 8.44E-04; 3.03E-03]; % Pa·s
Speedofsound = [1498; 1207; 1328; 1.04E+03; 1338]; % m/s
Molarvolume = [1.80E-05; 5.84E-05; 1.07E-04; 1.00E-04; 2.94E-04]; % m³/mol
Molecularvolume = [1.91E-29; 6.21E-29; 1.14E-28; 1.07E-28; 3.13E-28]; % m³

% Define row names correctly
RowNames = {
    'Solvent Name',
    'Molecular Weight [g/mol]', 
    'Molecular Mass [Kg]', 
    'Molecular Radius [m]', 
    'Density [kg/m³]', 
    'Viscosity [Pa·s]', 
    'Speed of Sound [m/s]', 
    'Molar Volume [m³/mol]', 
    'Molecular Volume [m³]'
};

% Create the table
solventproperties = table(SolventNamesRow, MW, solventmass, solventradius, Density, ...
    Viscosity, Speedofsound, Molarvolume, Molecularvolume, ...
    'VariableNames', RowNames, 'RowNames', SolventNames);

