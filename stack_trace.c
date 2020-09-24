#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>

/* gcc -g -rdynamic stack_trace.c -o stack_trace */
void stack_trace(){
  void *trace[16];
  char **messages = (char **)NULL;
  int i, trace_size = 0;

  trace_size = backtrace(trace, 16);
  messages = backtrace_symbols(trace, trace_size);
  printf("[stack trace]>>>\n");
  for (i=0; i < trace_size; i++)
    printf("%s\n", messages[i]);
  printf("<<<[stack trace]\n");
  free(messages);
}

void main(){
  stack_trace();
}
