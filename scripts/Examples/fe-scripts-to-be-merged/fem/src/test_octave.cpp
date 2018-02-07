#include <iostream>

#define OCTINTERP_API
#include <octave/octave.h>
#include <octave/oct.h>
#include <octave/parse.h>

using namespace std;

int main(int argc,char *argv[])
{
  int embedded;

  octave_main(argc,argv,embedded=1);

  FILE *file;
  file = fopen("/tmp/octave_empty.m","w");
  fprintf(file,"disp('running /tmp/octave_empty.m');\n");
  fclose(file);
  parse_and_execute ("/tmp/octave_empty.m");

  file = fopen("/tmp/print_x.m","w");
  fprintf(file,"x = 123.0;\ndisp(sprintf('x = %%e',x));\n");
  fclose(file);
  parse_and_execute("/tmp/print_x.m");

  file = fopen("/tmp/print_x2.m","w");
  fprintf(file,"disp(sprintf('x = %%e',x^2));\n");
  fclose(file);
  parse_and_execute("/tmp/print_x2.m");

  fprintf(stdout,"octave_debug = %d\n",octave_debug);
  fprintf(stdout,"line_editing = %d\n",line_editing);
  fprintf(stdout,"end_tokens_expected = %d\n",end_tokens_expected);
  // parse_and_execute (FILE *) does not work
  //FILE *istream;
  //istream = fopen("/tmp/print_x2.m","r");
  //parse_and_execute (istream);
  //fclose (istream);

  return embedded;
}
