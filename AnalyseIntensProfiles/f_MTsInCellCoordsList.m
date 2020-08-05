function [OutList] = f_MTsInCellCoordsList(Input)
OutList = cell(length(Input), 1);
for i_Cell = 1:length(Input)
    InCell = Input{i_Cell};
    for i = 1:length(InCell)
        CurrCell = InCell{i};
        OutList{i_Cell} = [OutList{i_Cell}; CurrCell(:, 1:2)];        
    end
end