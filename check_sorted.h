#define WALLCLOCK(time) do {                                   \
      unsigned long val;                                       \
      volatile unsigned int a, d;                              \
      __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d) : );   \
      val = ((unsigned long) a)|(((unsigned long)d)<<32);      \
      (time) = val / 3330000000.;                              \
    } while(0)


int check_sorted(int n, const mwIndex *ia, const mwIndex *ja)
{
    int i, j;

    for (i=0; i<n; i++)
    {
        for (j=ia[i]; j<ia[i+1]-1; j++)
        {
            if (ja[j] >= ja[j+1])
            {
                printf("%d not sorted %d %d\n", i, ja[j], ja[j+1]);
                return -1;
            }
        }
    }

    return 0;
}


