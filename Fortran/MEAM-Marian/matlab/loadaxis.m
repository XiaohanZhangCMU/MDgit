%computing stress tensor for loading axis

a = [ 1 -2  1 ]';
b = [ 1  0 -1 ]'; b2 = [ -1 1 0 ]'; b3 = [ 0 -1 1]';
c = [ 1  1  1 ]';

a = a/norm(a);
b = b/norm(b); b2 = b2/norm(b2); b3 = b3/norm(b3);
c = c/norm(c);

H = [a, b, c];
%H
disp(sprintf('det(H)=%f',det(H)));

%uniaxial compression direction
%axcub = [-2 9 20]';
%axcub = [ 2 9 20]';
%axcub = [-9 2 20]';
axcub = [-9 -2 20]';

axcub = axcub/norm(axcub);

burgcub = c;
normcub = b;
normcub2 = b2;
normcub3 = b3;


%stress in cubic axis
strcub = axcub*axcub';
schmidcub = burgcub'*strcub*normcub;
schmidcub2 = burgcub'*strcub*normcub2;
schmidcub3 = burgcub'*strcub*normcub3;

%stress in lab axis
strlab = H'*strcub*H;
burglab = H'*burgcub;
normlab = H'*normcub;
normlab2 = H'*normcub2;
normlab3 = H'*normcub3;

schmidlab = burglab'*strlab*normlab;
schmidlab2 = burglab'*strlab*normlab2;
schmidlab3 = burglab'*strlab*normlab3;


strcub,strlab

[schmidcub, schmidlab; schmidcub2, schmidlab2; schmidcub3, schmidlab3]