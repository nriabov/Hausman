function [ x ] = PrintFig( filo, fig )
%PRINTFIG: Takes in a figure handle and an output name and prints a PDF
%version of the figure in question

try
    print(fig,'-depsc2',sprintf('%s.eps',filo));
catch
    print(fig,'-depsc2',sprintf('%s.eps',filo));    
end

system(sprintf('ps2pdf -dEPSCrop %s.eps %s.pdf',filo,filo));
file_eps = sprintf('%s.eps',filo);
system(sprintf('del %s', file_eps));
x = 1;

end

