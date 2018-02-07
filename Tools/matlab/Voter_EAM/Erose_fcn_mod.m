function f = Erose_fcn_mod(a,a_c,epsilon)

%epsilon = Erose_fcn(a);
%epsilon_c = Erose_fcn(a_c);
f = (Erose_fcn((1-epsilon)^.5 * a) - epsilon)/(1-epsilon);
%f = (Erose_fcn((1-epsilon_c)^.5 * a) - epsilon_c)/(1-epsilon_c);
%f = (Erose_fcn((1-epsilon_c)^.5 * a) - epsilon)/(1-epsilon);

end