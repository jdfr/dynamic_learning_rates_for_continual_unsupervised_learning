FileName='Irregular.bmp';
NumSamples = 20000;
NumEpochs = 2;
NumSteps = NumSamples*NumEpochs;
NumRowsMap = 5;
NumColsMap = 5;

%Samples=rand(2,20000);
Samples = GenerateSamplesImg(FileName,NumSamples);

[Model] = TrainKohonenSOM(Samples,NumRowsMap,NumColsMap,NumSteps);

Handle = PlotSOM(Model);
h = plot(Samples(1,:),Samples(2,:),'*y');
uistack(h,'bottom');