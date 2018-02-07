#run with sw

MD++ setnolog
MD++ setoverwrite
MD++ dirname = "runs/test_octave"

puts "hello, world!!!"
puts "Print special characters: \$, \\, \[, \""

OCTAVE_Run /tmp/print_x.m 
OCTAVE_Run /tmp/print_x2.m 

OCTAVE disp('hello from octave!');
OCTAVE format long
OCTAVE { x = 123.0; }
OCTAVE { disp(sprintf('x   = %e',x));   }
OCTAVE { disp(sprintf('x^2 = %e',x^2)); }
OCTAVE { disp(sprintf('pi  = %e',pi));  }

OCTAVE "h=rand(3,3); h, inv(h)"
OCTAVE for i=1:4, disp(sprintf('i=%d',i))\; end 

MD++ {
  latticestructure = diamond-cubic
  latticeconst = 5.4309529817532409 #(A) for Si
  latticesize = [ 1  0  0 3    #(x)
                  0  1  1 3    #(y)
                  0  0  1 3  ] #(z)
  makecrystal finalcnfile = perf.cn writecn
}
MD++ eval

set x 0.199
OCTAVE x = $x
OCTAVE disp(sprintf('x = %f',x))\;

OCTAVE { h = zeros(3,3); }
for { set i 1 } { $i <=3 } { incr i } { 
  for { set j 1 } { $j <=3 } { incr j } { 
    OCTAVE h($i,$j) = [MD++_Get H( [expr ($i-1)*3+($j-1)] ) ] \;
  }
}
OCTAVE h
OCTAVE hinv = inv(h)
OCTAVE V = det(h)

set V [OCTAVE_Get V]
puts "V = $V"

puts "bye."
#OCTAVE exit

