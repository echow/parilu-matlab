%synchronous
%level 0, 1, 2
s=[
 0   418   418   418
 1   365   330   330
 2   325   265   250
 3   310   238   216
 4   299   222   191
 5   294   213   176
 6   290   206   165
 7   249   202   158
 8   286   199   152
 9   283   197   148
10   280   196   145
11   276   195   143
12   274   195   141
13   273   194   139
14   272   194   138
15   272   194   138
16   271   194   137
17   271   194   137
18   271   194   136
19   271   194   136
20   271   194   136
];
clf
plot(s(:,1),s(:,2));hold on
plot(s(:,1),s(:,3));
plot(s(:,1),s(:,4));
legend('level 0', 'level 1', 'level 2');
xlabel('Number of sweeps');
ylabel('PCG iterations');
xlim([0 15])
grid on
fig=gcf;fig.PaperUnits='inches';fig.PaperPosition=[0 0 3*1.618 3];

%asynchronous
% ordering and scheduling?
%level 0,1,2
a0=[
 0   418   418   418   418   418
 1   317   315   267   315   319
 2   272   275   266   268   271
 3   269   269   268   269   269
 4   270   271   271   271   271
 5   271   271   271   271   271
 6   271   271   271   271   271
];
a=a0;
clf;
plot(a(:,1),a(:,2));hold on
plot(a(:,1),a(:,3));
plot(a(:,1),a(:,4));
plot(a(:,1),a(:,5));
plot(a(:,1),a(:,6));
xlabel('Number of sweeps');
ylabel('PCG iterations');
grid on
fig=gcf;fig.PaperUnits='inches';fig.PaperPosition=[0 0 3*1.618/15*6 3];

a1=[
 0   418   418   418   418   418
 1   243   243   240   242   239
 2   199   199   199   198   199
 3   194   194   194   194   194
 4   194   194   194   194   194
 5   193   194   194   193   193
 6   193   193   193   193   193
];
a=a1;
clf;
plot(a(:,1),a(:,2));hold on
plot(a(:,1),a(:,3));
plot(a(:,1),a(:,4));
plot(a(:,1),a(:,5));
plot(a(:,1),a(:,6));
xlabel('Number of sweeps');
ylabel('PCG iterations');
grid on
fig=gcf;fig.PaperUnits='inches';fig.PaperPosition=[0 0 3*1.618 3];

a2=[
 0   418   418   418   418   418
 1   193   206   204   208   203
 2   146   146   145   145   145
 3   137   138   137   137   136
 4   135   135   135   135   135
 5   135   135   135   135   135
 6   135   135   135   135   135
];
a=a2;
clf;
plot(a(:,1),a(:,2));hold on
plot(a(:,1),a(:,3));
plot(a(:,1),a(:,4));
plot(a(:,1),a(:,5));
plot(a(:,1),a(:,6));
xlabel('Number of sweeps');
ylabel('PCG iterations');
grid on
fig=gcf;fig.PaperUnits='inches';fig.PaperPosition=[0 0 3*1.618/15*6 3];

%    level 0                level 1                level 2
a=[
  0  404  1.7e+04  41.1350  404  2.3e+04  41.1350  404  2.3e+04  41.1350
  1  318  3.8e+03  32.7491  256  5.7e+03  18.7110  206  7.0e+03  17.3239
  2  301  9.7e+02  32.1707  207  8.6e+02  12.4703  158  1.5e+03   6.7618
  3  298  1.7e+02  32.1117  193  1.8e+02  12.3845  132  4.8e+02   5.8985
  4  297  2.8e+01  32.1524  187  4.6e+01  12.4139  127  1.6e+02   5.8555
  5  297  4.4e+00  32.1613  186  1.4e+01  12.4230  126  6.5e+01   5.8706
%-1  297  0        32.1629  185  0        12.4272  126  0         5.8894
];

clf;
plot(a(:,1),a(:,2));hold on
plot(a(:,1),a(:,5));
plot(a(:,1),a(:,8));
xlabel('Number of sweeps');
%ylabel('PCG iterations');
grid on
fig=gcf;fig.PaperUnits='inches';fig.PaperPosition=[0 0 3*1.618/15*6 3];
