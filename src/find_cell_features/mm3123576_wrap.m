fprintf('-----MATLAB-BEGIN-----\n');
try
	rv = evalc('mm3123576');
	fprintf('-----SUCCESS\n%s', rv);
catch
	fprintf('-----ERROR\n%s\n%s', lasterr);
end
quit;
