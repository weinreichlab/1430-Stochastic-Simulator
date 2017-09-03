%this function creates the popup dialog that allows us to edit the mutation
%matrix
function new = InitialCountsDialog(current, numLoci, total, capacity)

    %Start Script
    d = dialog('Position',[100 100 400 400],'Name','Initial Counts');
    new = [];

    
    
    btn1 = uicontrol('Parent',d,...
               'Position',[150 25 100 25],...
               'String','Cancel',...
               'Callback',@cancel);
           
    btn2 = uicontrol('Parent',d,...
           'Position',[150 75 100 25],...
           'String','Save',...
           'Callback',@save);  
       
    btn3 = uicontrol('Parent',d,...
           'Position',[150 125 100 25],...
           'String','Equal Counts',...
           'Callback',@makeUniform);  
       
    
    widths = cell(1,size(current,1));
    for x = 1:size(current,1)
        widths{x} = 100;
    end
       
    
    table = uitable(gcf,...
        'Parent',d,...
        'Data', current(:),...
        'Position',[0 150 400 200],...
        'ColumnEditable', logical(ones(1,size(current,1))),...
        'ColumnWidth', widths);
    table.Data = num2cell(table.Data);
    rowname = cell(1, 2^numLoci);
    colname = {'Count'};
    for i = 1:(2^numLoci)
        rowname{i} = sprintf('Type %s', dec2bin(i - 1, numLoci));
    end
    
    table.RowName =  rowname;
    table.ColumnName =  colname;
    uiwait; %prevents the function from returning until uiresume is called

    %Functions
    
    %closes the window. uiresume 
    function cancel(~,~)
        uiresume;
        delete(gcf);
    end

    %Saves the data in the matrix to the matrix variable, which is returned
    %to the GUI script upon the execution of uiresume
    function save(~,~)
       data = cell2mat(table.Data);
       if any(isnan(data)) || ~isnumeric(data) || any(data < 0)
           warndlg('ERROR: All entries must be non-negative numbers!');
           return;
       elseif sum(data) ~= total && capacity
           warndlg(sprintf('ERROR: The sum of the counts must be %d!', total));
           return;
       elseif sum(data) > total && ~capacity
           warndlg(sprintf('ERROR: The sum of the counts must be no greater than %d!', total));
           return;
       end
       new = data';
       uiresume;
       delete(gcf);
    end

    function makeUniform(~,~)
        if capacity
            num = floor(total/length(table.Data));
            for j = 2:length(table.Data)
                table.Data{j} = num;
            end
            table.Data{1} = total - num*(length(table.Data) - 1);
        else
            for j = 1:length(table.Data)
                table.Data{j} = 1;
            end
        end
    end

    

end

