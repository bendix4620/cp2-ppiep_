#include "cfg.h"
#include <stdio.h>

#define CONFIG "../../configs_tests/cfg"

int readd(char *name)
{
    printf("value %s = %f\n", name, cfgGetd(name));
    return 0;
}
int readi(char *name)
{
    printf("value %s = %i\n", name, cfgGeti(name));
    return 0;
}

int main(void)
{
    cfgAlloc(CONFIG);
    printf("allocated\n");
    readd("test1");
    readd("test3");
    readd("test2");

    readi("test4");
    readi("test6");
    readi("test5");

    readd("test8");

    cfgFree();
    return 0;
}
