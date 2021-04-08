function observed = data_selection_three_groups(data, start_date, end_date)

% convert data type for date entries from string into datetime 
start_date = datetime(start_date);
end_date = datetime(end_date);

% index table by rows within date range
index = (data.Dates >= start_date) & (data.Dates <= end_date);

% subset table according to index
observed = data{index,{'Group1', 'Group2', 'Group3'}};

end