
#ifndef WIN_LINUX_AUX_H
#define WIN_LINUX_AUX_H 1

extern "C" void sleepMilli(int *milli);
extern "C" void pid(int *i);
extern "C" void hostname(char **h, int *i);

#ifdef WIN32 
#define SINCOS(phi, e0, e1) e0 = sin(phi); e1 = cos(phi);
#else
#define SINCOS(phi, e0, e1) sincos(phi, &(e0), &(e1));
#endif

#endif /* WIN_LINUX_AUX_H */


