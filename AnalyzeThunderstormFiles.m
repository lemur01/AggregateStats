%% Parameter
% Distance threshold
distTh = 90 % nm

% Synapse filters
LowerThSynSigma_nm_ = 80;
UpperThSynSigma_nm_ = 800;
LowerThSynIntensity_photon_ = 200;

input_folder =uigetdir() 

%% Analysis based on Thunderstorm output
mutant  = [];
exp = [];
ratio75 = []
at180 =[]
syn = []

folder = input_folder;
files=dir(fullfile(folder,'*0005.csv'));

fid = fopen(fullfile(folder, 'Summary.csv'), 'wt');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s, %s\n', ...
    'FWHMMin',  'FWHMMax', 'Eccentricity', 'Area', 'Shape','Amplitude', 'IntegratedInt','Theta' ,'Bg' , 'xAgg', 'yAgg','Vol','rmse', 'adjR', 'eFlag', 'CloseToSyn','DistToSyn','NrAgg', 'Filename');  % header
fclose(fid);

w = 5
for k = 1:length(files)
    smin = []; 
    smax = []; 
    ampl = []; 
    theta = []; 
    bg = []; 
    c1 = []; 
    c2 = []; 
    xAgg = [];yAgg = []; rmse = []; eflag = []; adj = [];
    yy = [];
    ch5 = readtable(fullfile(files(k).folder, files(k).name));
    Im = double(imread(strrep(fullfile(files(k).folder, files(k).name),'csv', 'tif')));

    xS = round(ch5.x_nm_/20); 
    yS = round(ch5.y_nm_/20); 
    ii = find((xS>w)&(yS>w)&(xS<=size(Im,2)-w)&(yS<=size(Im,1)-w));
    
    jj = setdiff(1:size(xS,1),ii);   
    fitS = arrayfun(@(i)( reshape(Im(yS(i)-w:yS(i)+w, xS(i)-w:xS(i)+w), [], 1 )), ii, 'UniformOutput',false);
   
    data = cell2mat(fitS');
    data = single(data)-2^16/2;
%   fit options
%   Lower,StartPoint and Upper. Each of the fields is an array of 7 values : 
%             offset
%             amplitude of the 2D gaussian
%             centroid X of the 2D gaussian
%             centroid Y of the 2D gaussian
%             angle of rotation for the 2D gaussian
%             width X of the 2D gaussian
%             width Y of the 2D gaussian

    fop.Lower = [0 0 w-1 w-1 0 0.2  0.2];
    fop.StartPoint = [0 max(data(:,1)) w+1 w+1 0 2/2.355  2/2.355];
    fop.Upper = [5000 5*max(data(:,1)) w+4 w+4 180 20 20 ];

    for ix = 1:size(data,2)
        ix
        [fitresult, gof, fout] = Gauss2DRotFit(reshape(data(:,ix), 11,11),fop)
        smin(ix) = min(fitresult.w1,fitresult.w2); 
        smax(ix) = max(fitresult.w1,fitresult.w2); 
        ampl(ix)=fitresult.b;
        theta(ix)=fitresult.t1;
        bg(ix)=fitresult.a;
        xAgg(ix)=fitresult.c1-w+xS(ix);
        yAgg(ix)=fitresult.c2-w+yS(ix);
        rmse(ix) = gof.rmse;
        adj(ix)= gof.adjrsquare;
        eflag(ix) = fout.exitflag;
    end

    % synapses left unmodified
    ch7 = readtable(fullfile(files(k).folder, strrep(files(k).name, '0005', '0007')));
    id = find((ch7.sigma_nm_<LowerThSynSigma_nm_)|(ch7.sigma_nm_>UpperThSynSigma_nm_)|(ch7.intensity_photon_<LowerThSynIntensity_photon_))
    ch7(id,:)=[];
    
    
    syn = [syn size(ch7.x_nm_,1)]
    ratio75= [ratio75 size(ch7.x_nm_,1)/size(ch5.x_nm_,1)];
    %nearest neighbor in ch5 for each query point in ch7
    [IdxA, D57] = knnsearch([ch5.x_nm_(ii), ch5.y_nm_(ii)], [ch7.x_nm_, ch7.y_nm_], 'K',2);
    %nearest neighbor in ch7 for each query point in ch5
    [IdxS, D75] = knnsearch([ch7.x_nm_, ch7.y_nm_],[ch5.x_nm_(ii), ch5.y_nm_(ii)], 'K',2);

    selAgg = find(D75(:,1)<= distTh); % select only aggregates closer that distTh centre to centre
    [GC,GR] = groupcounts(IdxS(selAgg,1)); 
    syn0 = size(ch7.x_nm_,1)-length(GR)
    newTabCol = zeros(height(ch7), 1);
    ch7.("AggregateNr") = newTabCol;
    ch7.("Area") = pi*ch7.sigma_nm_.^2;
    if isempty(selAgg)
        GR = 1;
        GC = 0;
        counts = syn0;
        centers = 0;
    else
        ch7.AggregateNr(GR) = GC;
        [counts centers] = hist(ch7.AggregateNr, 0:max(GC));
    end
    
    figure; bar(centers, counts); title('Aggregate distribution (per synapse)');
    xlabel('Nr aggregates')
    saveas(gcf,fullfile(files(k).folder, strrep(files(k).name, '0005.csv', '0007AggNr.png')))
    writetable(ch7,fullfile(files(k).folder, strrep(files(k).name, '0005', '0007Full')));

    Vol = 2*pi*ampl.*smin.*smax*400;
    Area = pi*smin.*smax*400;
    Shape = smin.*smax*400*2./(pi*(400*smin.^2+400*smax.^2));
    T = table(smin'*2.355*20,smax'*2.355*20,smin'./smax', Area', Shape', ampl', ampl'.*(2*pi).*smin'.*smax',theta',bg',xAgg'*20,yAgg'*20,Vol',rmse',adj',eflag',(D75(:,1)<= distTh),D75(:,1),...
        'VariableNames', {'FWHMMin',  'FWHMMax', 'Eccentricity', 'Area', 'Shape','Amplitude', 'IntegratedInt', 'Theta' ,'Bg' , 'xAgg', 'yAgg','Vol','rmse', 'adjR', 'eFlag', 'CloseToSyn','DistToSyn'});
    %writetable(T,fullfile(files(k).folder, strrep(files(k).name, '0005', '0005Full')), 'WriteMode','append');
    writetable(T,fullfile(files(k).folder, '0005FullResult.csv'), 'WriteMode','append');
    if size(table2array(T(T.CloseToSyn==1,:)),1) == 1
        summary = (T(T.CloseToSyn==1,:));
        
    else
        summary = mean(T(T.CloseToSyn==1,:));
    end
    
    aux = find(T.CloseToSyn==1);
    summary.("NrAgg") = size(aux,1);
    summary.("Filename") = files(k).name;
    writetable(summary,fullfile(files(k).folder, 'Summary.csv'),"WriteMode","append"); 
    
    if size(table2array(T(T.CloseToSyn==0,:)),1) == 1
        summary = (table2array(T(T.CloseToSyn==0,:)));
    else
        summary = mean((T(T.CloseToSyn==0,:)));
    end
    
    aux = find(T.CloseToSyn==0);
    summary.("NrAgg") = size(aux,1);
    summary.("Filename") = files(k).name;
    writetable(summary,fullfile(files(k).folder, 'Summary.csv'),"WriteMode","append");%,"AutoFitWidth",false);

    writetable(table(centers', counts','VariableNames', {'Values',  'Counts'} ),fullfile(files(k).folder, 'DistributionPerSynapse.csv'), 'WriteMode','append');

end
