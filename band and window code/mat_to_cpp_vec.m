function outstr = mat_to_cpp_vec(vec)
%takes a vector in matlab and converts it into a string of a c++ formatted 
%vector to then copy and paste into cpp code
%vec must contain numbers

outstr = '{';

for i = 1:length(vec)
    if i == length(vec)
        outstr = strcat(outstr, num2str(vec(i)), '}');
    else    
        outstr = strcat(outstr, num2str(vec(i)), ',');
    end
end

end