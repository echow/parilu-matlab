#include <stdlib.h>
#include <stdio.h>

void main()
{
    int n = 203841;

    int i;
    srand(123);

    for (i=0; i<n; i++)
    {
       printf("%.15f\n", rand()/(double)RAND_MAX - .5);
    }
}
