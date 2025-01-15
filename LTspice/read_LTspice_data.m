function tbl = read_LTspice_data(filename)

tbl0 = readtable(filename, 'VariableNamingRule', 'preserve');
tbl = table('Size',[size(tbl0.("Freq."),1) 1],'VariableTypes',{'double'}, 'VariableNames',{'freq'});

tbl.freq = tbl0.("Freq.");

counter = 1;
for name_cell = tbl0.Properties.VariableNames
    name = convertCharsToStrings(name_cell{1});
    if name ~= "Freq."
        data_cells = cellfun(@(x)sscanf(x, '(%fdB,%f'), tbl0.(name), 'Unif',0);
        data = cell2mat(data_cells.').';
        complex_data = 10.^(data(:, 1)/20) .* exp(1j*data(:, 2)*pi/180);
        tbl = addvars(tbl, complex_data, 'NewVariableNames', sprintf("G%u", counter));
        counter = counter + 1;
    end
end

end