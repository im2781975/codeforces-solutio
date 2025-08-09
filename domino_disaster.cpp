#include<iostream>
#include<string>
https://codeforces.com/problemset/problem/1367/A
using namespace std;

int main(){
    int t;
    int n;
    string s;

    cin>>t;
    for(int i = 0; i < t; i++){
      cin>>n;
      cin>>s;
      for(int j = 0; j < n; j++){
        if(s[j] == 'D'){

          s[j] = 'U';
          continue;
        }
        else if(s[j] == 'U'){

          s[j] = 'D';
          continue;

        }
        else{
          continue;
        }
      }
      cout<<s<<endl;

    }
    return 0;
}
