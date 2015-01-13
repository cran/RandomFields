
#include <math.h>
#include <unistd.h>
 
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include "RF.h"
#include <Rdefines.h>
//#include <curses.h>


/*

// http://stackoverflow.com/questions/421860/capture-characters-from-standard-input-without-waiting-for-enter-to-be-pressed

#include <sys/filio.h>
int kbhit()
{
 int i;
 ioctl(0, FIONREAD, &i);
 return i; // return a count of chars available to read 
}
main()
{
 int i = 0;
 intc='';
 system("stty raw -echo");
 //printf("enter 'q' to quit \n");
 for (;c!='q';i++) {
    if (kbhit()) {
        c=getchar();
       //printf("\n got %c, on iteration %d",c, i);
    }
}
 system("stty cooked echo");
}
*/


SEXP getChar() {
  error("does not work");
#ifdef WIN32
  error("input limitations on windows");
#endif
#define maxGetChar 255
  //typedef char intchar[sizeof(int) / sizeof(char)];
  //typedef char intchar[sizeof(int) / sizeof(char)];
  SEXP str;
  int //g,
    i = -1;
  char // c, *t = NULL,
    *s = NULL
    ;
  s = (char*) MALLOC(sizeof(char) * maxGetChar);
  // initscr();
  //  fflush(stdin);
  //  nocbreak();
  while (++i < maxGetChar) {
    // g = getchar();       
    //s[i] = ((intchar*) &g)[0][0];
    // g = scanf("%c\n", s); s[1] = '\0'; break;
    //t = fgets(s, 2, stdin); break;
    //    s[i] = getch();
    //fflush(stdin);
    if (false) {
      s[i+1] = '\0';
      // printf("%d i=%d  '%c' '%c' '%c' '%c' '%c'\n", g, i,  s[i],
      // ((intchar*) &g)[0][0],
      //  ((intchar*) &g)[0][1],
      // ((intchar*) &g)[0][2],
      // ((intchar*) &g)[0][3]
      // );
    }
    if (s[i] == '\n') {
      s[i] = '\0';
      break;
    }
  }
  //endwin();
//printf(">%s<\n", s);
  PROTECT(str=allocVector(STRSXP, 1));
  SET_STRING_ELT(str, 0, mkChar(s));  
  UNPROTECT(1);
  free(s);
  return str;
}
