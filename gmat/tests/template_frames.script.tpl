%----------------------------------------
%--   GMAT validate frame conversions  --
%----------------------------------------

% ABOUT
% -----
%
% The following script validates a particular
% frame conversion between a bodyICRS frame
% and its associated bobyFixed one.

%----------------------------------------
%---------- Coordinate Systems
%----------------------------------------

% Define the bodyICRS frame
Create CoordinateSystem bodyICRF;
GMAT bodyICRF.Origin = {{ body }};
GMAT bodyICRF.Axes = ICRF;

% Define the bodyFixed frame
Create CoordinateSystem bodyFixed;
GMAT bodyFixed.Origin = {{ body }};
GMAT bodyFixed.Axes = BodyFixed;

%----------------------------------------
%---------- Spacecraft
%----------------------------------------

% Build spacecraft's orbit from vectors
Create Spacecraft ss_0;
GMAT ss_0.DateFormat = TTGregorian;
GMAT ss_0.CoordinateSystem = bodyICRF;
GMAT ss_0.DisplayStateType = Cartesian;
GMAT ss_0.X = {{ rx }}; 
GMAT ss_0.Y = {{ ry }};
GMAT ss_0.Z = {{ rz }};
GMAT ss_0.VX = {{ vx }};
GMAT ss_0.VY = {{ vy }};
GMAT ss_0.VZ = {{ vz }};

%----------------------------------------
%---------- ForceModel
%----------------------------------------

% Build a two-body two-point force model
Create ForceModel PointMass;
GMAT PointMass.CentralBody = {{ body }};
GMAT PointMass.PointMasses = { {{ body }} };
GMAT PointMass.Drag = None;
GMAT PointMass.SRP = Off;
GMAT PointMass.RelativisticCorrection = Off;
GMAT PointMass.ErrorControl = RSSStep;


%----------------------------------------
%---------- Propagators
%----------------------------------------

% Create propgator with point mass only
Create Propagator bodyPointMass;
GMAT bodyPointMass.FM = PointMass;
GMAT bodyPointMass.Type = PrinceDormand78;
GMAT bodyPointMass.InitialStepSize = 0.5;
GMAT bodyPointMass.Accuracy = 1e-12;
GMAT bodyPointMass.MinStep = 0;
GMAT bodyPointMass.MaxStep = 86400;
GMAT bodyPointMass.MaxStepAttempts = 500;
GMAT bodyPointMass.StopIfAccuracyIsViolated = false;

%----------------------------------------
%---------- Solvers
%----------------------------------------

% Default GMAT solver
Create DifferentialCorrector DC;
GMAT DC.ShowProgress = true;
GMAT DC.ReportStyle = Normal;
GMAT DC.ReportFile = 'DifferentialCorrectorDC.data';
GMAT DC.MaximumIterations = 25;
GMAT DC.DerivativeMethod = ForwardDifference;
GMAT DC.Algorithm = NewtonRaphson;

%----------------------------------------
%---------- Plots/Reports
%----------------------------------------

Create OrbitView GLPlot;
GMAT GLPlot.Add = {ss_0, {{ body }}};
GMAT GLPlot.CoordinateSystem = bodyICRF;
GMAT GLPlot.ViewPointReference = {{ body }};
GMAT GLPlot.ViewDirection = {{ body }};

Create ReportFile report;
GMAT report.SolverIterations = Current;
GMAT report.Filename = 'ReportData.txt';
GMAT report.Precision = 7;
GMAT report.Add = {ss_0.TTGregorian, ss_0.ElapsedSecs, ss_0.bodyICRF.X, ss_0.bodyICRF.Y, ss_0.bodyICRF.Z, ss_0.bodyFixed.X, ss_0.bodyFixed.Y, ss_0.bodyFixed.Z};
GMAT report.WriteHeaders = true;
GMAT report.LeftJustify = Off;
GMAT report.ZeroFill = On;
GMAT report.FixedWidth = true;
GMAT report.Delimiter = ' ';
GMAT report.ColumnWidth = 10;
GMAT report.WriteReport = true;


%----------------------------------------
%---------- Mission Sequence
%----------------------------------------
BeginMissionSequence;

Propagate 'Propagate 1 hour' bodyPointMass(ss_0) {ss_0.ElapsedSecs = 3600};
