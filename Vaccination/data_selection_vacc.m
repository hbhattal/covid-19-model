function observed = data_selection_vacc(data, start_date)

first = datenum(start_date)-datenum('Jan 26, 2020')+2;

observed = data(first,:);

end