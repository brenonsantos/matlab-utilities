% Breno Niehues dos Santos - 2023/11/09


function export2pdf(filename)
% exports the latest figure to pdf
    clear figure_property;
    figure_property.units = 'inches';
    figure_property.format = 'pdf';
    figure_property.Preview= 'none';
    figure_property.Width= '5'; % Figure width on canvas
    figure_property.Height= '4'; % Figure height on canvas
    figure_property.Units= 'inches';
    figure_property.Color= 'rgb';
    figure_property.Background= 'w';
    figure_property.FixedfontSize= '12';
    figure_property.ScaledfontSize= 'auto';
    figure_property.FontMode= 'scaled';
    figure_property.FontSizeMin= 'auto';
    figure_property.FixedLineWidth= '1';
    figure_property.ScaledLineWidth= 'auto';
    figure_property.LineMode= 'none';
    figure_property.LineWidthMin= '1';
    figure_property.FontName= 'Helvetica';% Might want to change this to something that is available
    figure_property.FontWeight= 'auto';
    figure_property.FontAngle= 'auto';
    figure_property.FontEncoding= 'latin1';
    figure_property.PSLevel= '3';
    figure_property.Renderer= 'painters';
    figure_property.Resolution= '600';
    figure_property.LineStyleMap= 'none';
    figure_property.ApplyStyle= '0';
    figure_property.Bounds= 'tight';
    figure_property.LockAxes= 'off';
    figure_property.LockAxesTicks= 'off';
    figure_property.ShowUI= 'off';
    figure_property.SeparateText= 'off';
    chosen_figure=gcf;
    set(chosen_figure,'PaperUnits','inches');
    set(chosen_figure,'PaperPositionMode','auto');
    set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
    set(chosen_figure,'Units','inches');
    
    filename = strcat(filename, '.pdf');
    hgexport(gcf, filename, figure_property); %Set desired file name
end

