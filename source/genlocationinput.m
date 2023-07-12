function genlocationinput(path, nLine, dims)
if nargin < 3, dims = [20, 20]; end;
if nargin < 2, nLine = 70589; end;
if nargin < 1, path = '../input/location.input'; end;
nRow = dims(1);
nCol = dims(2);
n = 1;
x = 1;
y = 1;
f = fopen(path, 'w');
for i = 1 : nLine
    fprintf(f, '%d\t%d\t%d\n', n, x, y);
    y = y + 1;
    if y > nCol, 
        y = 1;
        x = x + 1;
        if x > nRow,
            x = 1;
            n = n + 1;
        end;
    end;
end;
fclose(f);
    
