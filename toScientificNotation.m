function sciStr = toScientificNotation(num)
    % Convert the number to scientific notation (mantissa and exponent)
    [exponent] = log10(abs(num));
    
    % Round mantissa to have only three significant digits
    coefficient = num / 10^floor(exponent);
    
    % Construct the final string in the format X.XXdY
    sciStr = sprintf('%.4fd%d', coefficient, floor(exponent));
end
