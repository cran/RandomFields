
#ifndef WIN_LINUX_AUX_H
#define WIN_LINUX_AUX_H 1

extern "C" void sleepMilli(int *milli);
extern "C" void sleepMicro(int *milli);
extern "C" void pid(int *i);
extern "C" void hostname(char **h, int *i);


#endif /* WIN_LINUX_AUX_H */


