% Conversion script
% This function converts all the .fig files to .pdf

addpath 'Y:\McCloskey2014-Sims\Hausman\Scripts'

beta1_vec = 0.1:0.1:0.9;
%dirname = 'Y:\McCloskey2014-Sims\Hausman\beta1search\';
dirname = 'Y:\McCloskey2014-Sims\Hausman\beta1search-twosided\';


for i = 1:length(beta1_vec)
%for i = 1:1
    cur_beta1 = beta1_vec(i);
    curdir = ['beta1-',num2str(cur_beta1)];
    %curdir = 'pdf_new'
    cd([dirname, curdir]);
    mkdir('pdf_new');
    %cd('pdf_new');
    %mkdir([curdir, '-pdf'])

    
    %for j = 1:1
    for j = 1:8
        % Set figure name
        %name = ['PowerCurves-twosided-h2-', num2str(j)];
        name = ['PowerCurves-h2-', num2str(j)];

        % Load the figure
        %load([name, '.fig'], '-mat')
        openfig([name, '.fig'], 'reuse', 'invisible');
        
        % Output figure as .pdf
        %cd([curdir, '-pdf'])
        cd('pdf_new');
        PrintFig(name, gcf);
        clf
        close all
        cd ..
    end
end