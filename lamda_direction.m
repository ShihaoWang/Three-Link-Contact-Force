function lamda = lamda_direction(Theta)

i = [ 1 0 0 ]';

j = [ 0 1 0 ]';

lamda = cos(Theta) * i + sin(Theta) * j;

end