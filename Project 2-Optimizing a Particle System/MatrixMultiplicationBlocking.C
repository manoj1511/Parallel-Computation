#include <iostream>
#include <stdlib.h>
#include <ctime>
using namespace std;
#define N 1000

int main(int argc, const char** argv)
{
  clock_t start_time,end_time;
  start_time = clock();
  int BLOCK_SIZE = 73;
  int n = BLOCK_SIZE * (N/BLOCK_SIZE);
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
  int sum=0;
  for(int i1=0;i1<n;i1+=BLOCK_SIZE)
  {
      for(int j1=0;j1<n;j1+=BLOCK_SIZE)
      {
          for(int k1=0;k1<n;k1+=BLOCK_SIZE)
          {
              for(int i=i1;i<i1+BLOCK_SIZE;i++)
              {
                  for(int j=j1;j<j1+BLOCK_SIZE;j++)
                  {
                      sum = c[i][j];
                      for(int k=k1;k<k1+BLOCK_SIZE;k++)
                      {
                          sum += a[i][k] * b[k][j];
                      }
                      c[i][j] = sum; 
                  }
                 
              }
          }
      }
         }
        cout << "\nthe array size is :"<<N<<"\n";
        cout << "the Block size is :"<<BLOCK_SIZE<<"\n";
        cout << "the final array has been created\n";
        cout << "time to see the performance \n";
        end_time = clock();
        cout<<"Time Taken is : "<<(end_time - start_time)/1000<<" ms"<<endl;
       /*cout << "\nthe final array has been created\n\n";
        cout << "time to see the results \n";
       for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                cout<<"\t"<<c[i][j];
            }
            cout<<"\n\n";
        }*/
          
  return 0;
}