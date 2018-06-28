#include <iostream>
#include <stdlib.h>
#include <ctime>
using namespace std;
#define N 2000

int main(int argc, const char** argv)
{
  clock_t start_time,end_time;
  start_time = clock();
  int n = N;
  static int a[N][N] ={0};
  static int b[N][N] ={0};
  static int c[N][N] ={0};
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            a[i][j] = rand()%100;
            b[i][j] = rand()%100;
        }
    }
 
  static int tmp[N][N] ={0};
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      tmp[i][j] = b[j][i];
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      for (int k = 0; k < N; ++k)
        c[i][j] += a[i][k] * tmp[j][k];
      
        cout << "\nthe array size is :"<<N<<"\n";
        cout << "the final array has been created\n";
        cout << "time to see the performance \n";
        end_time = clock();
        cout<<"Time Taken is : "<<(end_time - start_time)/1000<<" ms"<<endl;
      
          
  return 0;
}


 




