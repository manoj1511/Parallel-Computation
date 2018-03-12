#include <iostream>
#include <stdlib.h>
#include <ctime>
using namespace std;
#define N 100

int main(int argc, const char** argv)
{
  clock_t start_time,end_time;
  start_time = clock();
  int n = N;
  static double a[N][N] ={0};
  static double b[N][N] ={0};
  static double c[N][N] ={0};
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            a[i][j] = rand()%100;
            b[i][j] = rand()%100;
        }
    }
  int sum=0;
  
   
              for(int i=0;i<n;i++)
              {
                  for(int j=0;j<n;j++)
                  {
                      sum = c[i][j];
                      for(int k=0;k<n;k++)
                      {
                          sum += a[i][k] * b[k][j];
                      }
                      c[i][j] = sum; 
                  }
                 
              }
       
      
        cout << "\nthe array size is :"<<N<<"\n";
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