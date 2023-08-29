function fileDifference(file1Name, file2Name, outputFileName)


% Read the content of the files
file1Content = fileread(file1Name);
file2Content = fileread(file2Name);

% Split the content into lines
file1Lines = strsplit(file1Content, '\n');
file2Lines = strsplit(file2Content, '\n');

% Initialize a cell array to store the differences
differences = {};

% Find lines in file1 that are not present in file2
for lineNumber = 1:length(file1Lines)
    currentLine = file1Lines{lineNumber};
    if ~any(strcmp(currentLine, file2Lines))
        differences{end+1} = ['Line ' num2str(lineNumber) ': ' currentLine];
    end
end

% Write the differences to the output file
outputFile = fopen(outputFileName, 'w');
if isempty(differences)
    fprintf(outputFile, 'No differences found.');
else
    fprintf(outputFile, 'Lines in %s that are not in %s:\n', file1Name, file2Name);
    for i = 1:numel(differences)
        fprintf(outputFile, '%s\n', differences{i});
    end
end
fclose(outputFile);

disp(['Differences written to ' outputFileName]);
