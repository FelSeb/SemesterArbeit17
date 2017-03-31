function [ val, out2 ] = alternativelogpdf( x )


if(length(x) == 2)
    val(1,1) = log(normpdf(x(1),1,1));
    val(2,1) = log(normpdf(x(2),1,1));
    out2{1} = ['hi1  ', num2str(x(1))];
    out2{2} = ['hi2  ', num2str(x(2))];
else
    val = log(normpdf(x,1,1));
    out2 = 'hello';
end


end

