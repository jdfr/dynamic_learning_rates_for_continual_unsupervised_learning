function [FormattedString]=PrettyNumbers(Average,StdDeviation,Significant)

if abs(Average)>=0.1
    FormattedString=sprintf('%2.2f (%2.2f)',Average,StdDeviation);
    FormattedString=strrep(FormattedString,'-0.','-.');
    if strcmp(FormattedString(1:2),'0.')==1
        FormattedString=FormattedString(2:end);
    end
    FormattedString=strrep(FormattedString,'(0.','(.');
else
    FormattedString=sprintf('%2.2e (%2.2e)',Average,StdDeviation);
end

if Significant
    FormattedString=sprintf('%s*',FormattedString);
end