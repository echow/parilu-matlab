#include <stdlib.h>
#include <stdio.h>

void main()
{
    int n = 203841; // set to whatever length you need

    int i;
    srand(0); // set seed to 0

    for (i=0; i<n; i++)
    {
       printf("%.15f\n", rand()/(double)RAND_MAX - .5);
    }
}
