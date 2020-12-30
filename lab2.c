#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#define N 4938921   //基因组一共70*70556=4938920长 + $

//后缀结构体
typedef struct suffix{
    int SA;  //该后缀的SA
    char BWT;   //该后缀的BWT索引值
    int index;    //该后缀的索引
    char* str;  //该后缀 具体是串开始的地址
 }suffix;//后缀

 typedef struct SA_interval{
     int up[10000];     //所有上界
     int down[10000];   //所有下界 一一对应
     int index;         //现在用到了那个下标了
 }SA_interval;  //不精确比对时的

//全局变量部分
char* T = NULL;   //文本串
char* reverse_T = NULL;   //反转后的文本串
suffix* suf = NULL;    //所有后缀
suffix* reverse_suf = NULL;    //反转后的后缀
char* BWT = NULL;   //BWT索引
char* reverse_BWT = NULL;   //反转后的BWT
char reads1[10][75];   //文件一的read
char reads2[10][75];   //文件二的read
int count[10] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};    //分别标识$ A区  C  G  T的开始位置和结束位置   前两个一定是0 0
int appear[4][4938921] = {0};
int reverse_appear[4][4938921] = {0};
int D[70];      //不精确匹配时的D数组  生命值
SA_interval inEXA_intervals;   //记录当前匹配的所有区间

//获得整个参考基因组
int get_T()
{
    // char T[4938920]; //基因组
    T = (char*)malloc(5000000);

    FILE *fp = NULL;
    int i = 0,j;
    char c;

    fp = fopen("NC_008253.fna", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
    if(fp == NULL)
    {
        printf("打开文件失败");
        return -1;;
    }

    //读头
    for (i = 0;i < 69;i++)
    {
        fscanf(fp,"%c",&c);
    }

    //读基因组
    i = 0;
    for (j = 0;j < 5009476;j++)
    {
        fscanf(fp,"%c",&c);
        if (c == '\n')
            continue;
        else if (c != EOF)
        {
            T[i] = c;
            i++;
        }
        else
            break;        
    }    
    T[i] = '$';
    i++;
    T[i] = '\0';

    fclose(fp);
    return 1;
}

//提取出每一个后缀
suffix* get_suf()
{
    int j;
    suffix* suf = (suffix*)malloc(N*sizeof(suffix));//一共N个后缀 算上$

    //第一个后缀的BWT指定是最后一个字符的位置
    suf[0].SA = 0;   
    suf[0].BWT = T[N-1];   //最后一个字符
    suf[0].str = T;     //第一个后缀的地址

    //从第二个字符开始分配SA值和BWT值
    for (j = 1;j < N;j++)
    {        
        suf[j].SA = j; //几后缀  
        suf[j].BWT = T[j-1];    //每个后缀的BWT是它对应的SA的前一个字符
        suf[j].index = j;
        suf[j].str = T + j;
    }

    return suf;
}

//比较两个串的大小
int compare_string(const void * a, const void * b)  
{    
    suffix* aa = (suffix*)a;
    suffix* bb = (suffix*)b;
    int i = 0;
    while (1)
    {
        if (aa->str[i] == '$' && bb->str[i] == '$')  //和自己比  返回0
                return 0;
        else if (aa->str[i] == '$')   //有$出现就是短的比较小
            return -1;
        else if(bb->str[i] == '$')
            return 1;
        else if (aa->str[i] == bb->str[i])  //如果两字符相等 比不出来 看下一位
        {
            i++;
            continue;
        }
        else //两字符不相等而且没有是$符的  能比出来
        {
            if(aa->str[i] > bb->str[i])
                return 1;
            else
                return -1;
        }
        
    }
}

void sort_suf(suffix* suf)
{    
    qsort(suf, N, sizeof(suffix), compare_string);
}


void get_BWT()
{
    int i;
    BWT = (char*)malloc(N*sizeof(char));    //申请BWT大小的数组

    //将整个文本串读进内存
    get_T();
    //得到后缀及对应的BWT
    suf = get_suf();
    sort_suf(suf);         //得到SA及正确顺序的BWT

    for (i = 0;i < N;i++)   //将得到的BWT从结构体数组中拿出来
    {
        BWT[i] = suf[i].BWT;
    }
    
    printf("0");
}

void get_count_and_O()
{
    int i,j;
    int count_$ = 0, count_A = 0, count_C = 0, count_G = 0, count_T = 0;

    for (i = 0;i < N;i++)
    {
        if (i == 0)     //第一个BWT是C
        {
            appear[0][i] = 0; 
            appear[1][i] = 1;
            appear[2][i] = 0; 
            appear[3][i] = 0;
        }
        else
        {
            if (suf[i].BWT == 'A')
            {
                appear[0][i] = appear[0][i-1] + 1;      //当该字符是A 出现次数就加一 
                appear[1][i] = appear[1][i-1];
                appear[2][i] = appear[2][i-1]; 
                appear[3][i] = appear[3][i-1];
            }                            
            else if (suf[i].BWT == 'C')
            {
                appear[0][i] = appear[0][i-1];      //当该字符是C 出现次数就加一 
                appear[1][i] = appear[1][i-1] + 1;
                appear[2][i] = appear[2][i-1]; 
                appear[3][i] = appear[3][i-1];
            }             
            else if (suf[i].BWT == 'G')
            {
                appear[0][i] = appear[0][i-1];      //当该字符是G 出现次数就加一 
                appear[1][i] = appear[1][i-1];
                appear[2][i] = appear[2][i-1] + 1; 
                appear[3][i] = appear[3][i-1];
            }           
            else if (suf[i].BWT == 'T')
            {
                appear[0][i] = appear[0][i-1];      //当该字符是T 出现次数就加一 
                appear[1][i] = appear[1][i-1];
                appear[2][i] = appear[2][i-1]; 
                appear[3][i] = appear[3][i-1] + 1;
            }
        }
   


        //统计每个区的开始和结束
        if (suf[i].str[0] == '$')
            continue;
        else if (suf[i].str[0] == 'A')
            count_A++;       
        else if (suf[i].str[0] == 'C')
            count_C++;         
        else if (suf[i].str[0] == 'G')
            count_G++;         
        else if (suf[i].str[0] == 'T')
            count_T++;   
    }


    count[0] = 0;//不能变就是在0位置开始，在0位置结束
    count[1] = 0;   //表示$

    //更新A的
    if (count_A == 0)  
    {
        count[2] = -1;
        count[3] = -1;
    }          
    else 
    {
        count[2] = 1;
        count[3] = count_A;
    }

    //更新C的
    if (count_C == 0)
    {
        count[4] = -1;
        count[5] = -1;
    }        
    else 
    {
        count[4] = count_A + 1;
        count[5] = count_C + count_A;

    }
    
    //更新G的
    if (count_G == 0)
    {
        count[6] = -1;
        count[7] = -1;
    }        
    else
    {
        count[6] = count_C + count_A + 1;
        count[7] = count_G + count_C + count_A;
    } 

    //更新T的
    if (count_T == 0)
    {
        count[8] = -1;
        count[9] = -1;
    }        
    else 
    {
        count[8] = count_G + count_C + count_A + 1;
        count[9] = count_T + count_G + count_C + count_A;
    }
}

void get_reads()
{
    FILE *fp1 = NULL;   //第一个read文件
    FILE *fp2 = NULL;   //第二个read文件    

    int i = 0,j,index = 0,len;
    char c;
    char buf1[1024];  /*缓冲区*/
    char buf2[1024];  /*缓冲区*/   


    fp1 = fopen("Ecoli_4x.fq1", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
    fp2 = fopen("Ecoli_4x.fq2", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式

    if(fp1 == NULL ||fp2 == NULL)
    {
        printf("打开文件失败");
        return;
    }

    //先将文件一和文件二的read读入内存  这里度读10行
    for (i = 0;i < 40;i++)
    {        
        fgets(buf1,1024,fp1);      //按行读取文件1
        fgets(buf2,1024,fp2);       //按行读取文件2
        len = strlen(buf1);         //得到本行的长度
        buf1[len-1] = '\0';  //去掉换行符
        buf2[len-1] = '\0';  //去掉换行符
        if (i % 4 != 1)     //观察文件可知，第二行是reads，每4行是一个reads
            continue;
        else
        {
            memcpy(reads1[index],buf1,sizeof(char)*71);     //将缓冲区的read写入数组
            memcpy(reads2[index],buf2,sizeof(char)*71);     //将缓冲区的read写入数组
            index++;    //索引加一
        }
    }
   


    fclose(fp1);
    fclose(fp2);
}

///////////////////////////////////////exact//////////////////////////////////////////////////////
void exact_write_file(int R_up, int R_down, int file_num, int num, int direction)
{
    FILE *fp1 = NULL;
    FILE *fp2 = NULL;

    if (file_num = 1)   //比对的是文件一   写入文件一对应的SAM文件
    {
        fp1 = fopen("Alignment1.sam", "a"); //第一个逗号前是文件位置。逗号之后是打开文件方式
        if(fp1 == NULL)
        {
            printf("打开文件失败");
            return ;
        }

        fprintf(fp1,"\n%d\t%d\t%d\t", num, R_up, R_down);
        printf("\n%d\t%d\t%d\t", num, R_up, R_down);
        // printf("\n%d\t%d\t%d\t", num, R_up, R_down);
    }
    else        //比对的是文件二   写入文件二对应的SAM文件
    {
        fp1 = fopen("Alignment2.sam", "a"); //第一个逗号前是文件位置。逗号之后是打开文件方式
        if(fp1 == NULL)
        {
            printf("打开文件失败");
            return ;
        }

        fprintf(fp1,"\n%d\t%d\t%d\t", num, R_up, R_down);
        printf("\n%d\t%d\t%d\t", num, R_up, R_down);
        // printf("\n%d\t%d\t%d\t", num, R_up, R_down);
        // printf("%d\n", R_down);
    }

    fclose(fp1);
    fclose(fp2);
}

int find_exact_SA_interval(int num, int file_num)
{
    //对于第file_num个文件所有reads中第num个read    找对应的区间
    //写入文件   返回是否找到合法区间 1为找到   -1  没找到

    int i,j;
    int R_down,R_up;    //区间上下界  
  

    char read[70];  //一个read70个字符 

    if (file_num == 1)
    {
        memcpy(read,reads1[num],sizeof(char)*70);    //得到要比对的文件一的第num条read 
    }
    else
    {
        memcpy(read,reads2[num],sizeof(char)*70);   //得到要比对的文件二的第num条read 
    }
    

    for (i = 69;i >= 0;i--)     //从右往左比
    {
        
        if (read[i] == 'A')  //如果该位置是A
        {
            if(i == 69)     //如果是右边第一个字符的话
            {
                R_up = count[2];        //上界直接等于count开始
                R_down = count[3];      //下界直接等于结束的位置
            }
            else
            {
                R_up = count[2]-1 + appear[0][R_up-1] + 1;       //如果不是右边第一个字符 就按照正常迭代着来
                R_down = count[2]-1 + appear[0][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }            
        }
        else if (read[i] == 'C')
        {
            if(i == 69)
            {
                R_up = count[4];
                R_down = count[5];
            }
            else
            {
                R_up = count[4]-1 + appear[1][R_up-1] + 1;
                R_down = count[4]-1 + appear[1][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }
        }
        else if (read[i] == 'G')
        {
            if(i == 69)
            {
                R_up = count[6];
                R_down = count[7];
            }
            else
            {
                R_up = count[6]-1 + appear[2][R_up-1] + 1;
                R_down = count[6]-1 + appear[2][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }
        }
        else if (read[i] =='T')
        {
            if(i == 69)
            {
                R_up = count[8];
                R_down = count[9];
            }
            else
            {
                R_up = count[8]-1 + appear[3][R_up-1] + 1;
                R_down = count[8]-1 + appear[3][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }
        }
    }   //循环结束后得到70前缀的区间

    exact_write_file(R_up, R_down, file_num, num, 1);   //写入文件   最后的1代表的是正向比对的结果
    return 1;       //成功找到对应的区间  返回1   
}

int reverse_find_exact_SA_interval(int num, int file_num)
{
    //得到反向互补  再从左到右比对
    //返回是否找到合法区间

    int i,j;
    int R_down,R_up;    //区间上下界

    char orig_read[70];   //原read
    char read[70];      //反向互补后的read

    if (file_num == 1)  //得到原read
    {
        memcpy(orig_read,reads1[num],sizeof(char)*70);
    }
    else
    {
        memcpy(orig_read,reads2[num],sizeof(char)*70);
    }
    

    for (j = 69;j >= 0;j--) //得到要比对的read的反向互补
    {
        if(orig_read[j] == 'A')
        {
            read[70-j-1] = 'T';
        }
        else if(orig_read[j] == 'C')
        {
            read[70-j-1] = 'G';
        }
        else if(orig_read[j] == 'G')
        {
            read[70-j-1] = 'C';
        }
        else if(orig_read[j] == 'T')
        {
            read[70-j-1] = 'A';
        }
    }    

    //开始比对
    for (i = 69;i >= 0;i--)     //从右往左比
    {            
        if (read[i] == 'A')  //如果该位置是A
        {
            if(i == 69)     //如果是右边第一个字符的话
            {
                R_up = count[2];        //上界直接等于count开始
                R_down = count[3];      //下界直接等于结束的位置
            }
            else
            {
                R_up = count[2]-1 + appear[0][R_up-1] + 1;       //如果不是右边第一个字符 就按照正常迭代着来
                R_down = count[2]-1 + appear[0][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }            
        }
        else if (read[i] == 'C')
        {
            if(i == 69)
            {
                R_up = count[4];
                R_down = count[5];
            }
            else
            {
                R_up = count[4]-1 + appear[1][R_up-1] + 1;
                R_down = count[4]-1 + appear[1][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }
        }
        else if (read[i] == 'G')
        {
            if(i == 69)
            {
                R_up = count[6];
                R_down = count[7];
            }
            else
            {
                R_up = count[6]-1 + appear[2][R_up-1] + 1;
                R_down = count[6]-1 + appear[2][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }
        }
        else if (read[i] =='T')
        {
            if(i == 69)
            {
                R_up = count[8];
                R_down = count[9];
            }
            else
            {
                R_up = count[8]-1 + appear[3][R_up-1] + 1;
                R_down = count[8]-1 + appear[3][R_down];
                if (R_up > R_down)      //如果找不到对应区间  说明不是基因组的子串
                    return -1;
            }
        }
    }

    exact_write_file(R_up, R_down, file_num, num, -1);  //将比对结果写入文件   最后的1代表的是正向比对的结果  
    return 1;   
}

void exact_aligment()
{
    int i,j,flag;

    get_reads();    //将要比对的read读入内存

    //先按文件一比对  如果找不到区间的话  就求反向互补比对
    for (i = 0;i < 10;i++)
    {
        flag = find_exact_SA_interval(i, 1);  //传入read序号   先正向找
        if (flag = -1)
        {
            flag = reverse_find_exact_SA_interval(i, 1);   //正向找不到  找反向互补
        }
    }

    //再比文件二  先正向找  如果不正确的话再反向互补找
    for (i = 0;i < 10;i++)
    {
        flag = find_exact_SA_interval(i, 2);  //传入read序号   先正向找
        if (flag = -1)
        {
            flag = reverse_find_exact_SA_interval(i, 2);   //正向找不到  找反向互补
        }
    }
}

///////////////////////////////////////inexact////////////////////////////////////////////////////
void get_reverse_T()
{
    //得到反转的基因组
    int i;
    reverse_T = (char*)malloc(N);   //申请内存

    for (i = N-2;i >= 0;i--)
    {
        reverse_T[N-i-2] = T[i];    //对应位置复制
    }
    reverse_T[N-1] = T[N-1];

}

suffix* get_reverse_suf()
{
    //得到反转基因组的所有的后缀

    int j;
    suffix* suf = (suffix*)malloc(N*sizeof(suffix));//一共N个后缀 算上$

    //第一个后缀的BWT指定是最后一个字符的位置
    suf[0].SA = 0;   
    suf[0].BWT = reverse_T[N-1];   //最后一个字符
    suf[0].str = reverse_T;     //第一个后缀的地址

    //从第二个字符开始分配SA值和BWT值
    for (j = 1;j < N;j++)
    {        
        suf[j].SA = j; //几后缀  
        suf[j].BWT = reverse_T[j-1];    //每个后缀的BWT是它对应的SA的前一个字符
        suf[j].index = j;
        suf[j].str = reverse_T + j;
    }

    return suf;
}

void get_reverse_BWT()
{
    //得到反转基因组的BWT
    int i;
    reverse_BWT = (char*)malloc(N*sizeof(char));    //申请BWT大小的数组

    //得到反转的文本
    get_reverse_T();
    //得到后缀及对应的BWT
    reverse_suf = get_reverse_suf();
    sort_suf(reverse_suf);         //得到SA及正确顺序的BWT

    for (i = 0;i < N;i++)   //将得到的BWT从结构体数组中拿出来
    {
        reverse_BWT[i] = reverse_suf[i].BWT;
    }    
    printf("0");
}

void get_reverse_O()
{
    //得到反转基因组的appear函数 即O数组
    int i,j;

    for (i = 0;i < N;i++)
    {
        if (i == 0)     //第一个BWT是C
        {
            reverse_appear[0][i] = 0; 
            reverse_appear[1][i] = 1;
            reverse_appear[2][i] = 0; 
            reverse_appear[3][i] = 0;
        }
        else
        {
            if (reverse_BWT[i] == 'A')
            {
                reverse_appear[0][i] = reverse_appear[0][i-1] + 1;      //当该字符是A 出现次数就加一 
                reverse_appear[1][i] = reverse_appear[1][i-1];
                reverse_appear[2][i] = reverse_appear[2][i-1]; 
                reverse_appear[3][i] = reverse_appear[3][i-1];
            }                            
            else if (reverse_BWT[i] == 'C')
            {
                reverse_appear[0][i] = reverse_appear[0][i-1];      //当该字符是C 出现次数就加一 
                reverse_appear[1][i] = reverse_appear[1][i-1] + 1;
                reverse_appear[2][i] = reverse_appear[2][i-1]; 
                reverse_appear[3][i] = reverse_appear[3][i-1];
            }             
            else if (reverse_BWT[i] == 'G')
            {
                reverse_appear[0][i] = reverse_appear[0][i-1];      //当该字符是G 出现次数就加一 
                reverse_appear[1][i] = reverse_appear[1][i-1];
                reverse_appear[2][i] = reverse_appear[2][i-1] + 1; 
                reverse_appear[3][i] = reverse_appear[3][i-1];
            }           
            else if (reverse_BWT[i] == 'T')
            {
                reverse_appear[0][i] = reverse_appear[0][i-1];      //当该字符是T 出现次数就加一 
                reverse_appear[1][i] = reverse_appear[1][i-1];
                reverse_appear[2][i] = reverse_appear[2][i-1]; 
                reverse_appear[3][i] = reverse_appear[3][i-1] + 1;
            }
        }
    }
}

void INexact_write_file(int R_up, int R_down, int file_num, int num, int direction)
{
    //写入文件
    FILE *fp1 = NULL;
    FILE *fp2 = NULL;

    if (file_num = 1)   //比对的是文件一   写入文件一对应的SAM文件
    {
        fp1 = fopen("Inexact_Alignment1.sam", "a"); //第一个逗号前是文件位置。逗号之后是打开文件方式 以追加的方式
        if(fp1 == NULL)
        {
            printf("打开文件失败");
            return ;
        }

        fprintf(fp1,"\n%d\t%d\t%d\t", num, R_up, R_down);
        // printf("%d\n", R_up);
        // printf("%d\n", R_down);
    }
    else        //比对的是文件二   写入文件二对应的SAM文件
    {
        fp1 = fopen("Inexact_Alignment2.sam", "a"); //第一个逗号前是文件位置。逗号之后是打开文件方式 以追加的方式
        if(fp1 == NULL)
        {
            printf("打开文件失败");
            return ;
        }

        fprintf(fp1,"\n%d\t%d\t%d\t", num, R_up, R_down);
        // printf("%d\n", R_up);
        // printf("%d\n", R_down);
    }

    fclose(fp1);
    fclose(fp2);
}

void calculate_D(char read[])
{
    //不精确匹配时  计算D数组
    int k,l,z,i;

    k = 1;
    l = N-1;
    z = 0;

    for (i = 0;i < 70;i++)      //从左往右比对
    {
        if (read[i] == 'A')
        {
            k = count[2]-1 + reverse_appear[0][k-1] + 1;
            l = count[2]-1 + reverse_appear[0][l];
        }        
        else if (read[i] == 'C')
        {
            k = count[4]-1 + reverse_appear[1][k-1] + 1;
            l = count[4]-1 + reverse_appear[1][l];
        }
        else if (read[i] == 'G')
        {
            k = count[6]-1 + reverse_appear[2][k-1] + 1;
            l = count[6]-1 + reverse_appear[2][l];
        }
        else if (read[i] == 'T') 
        {
            k = count[8]-1 + reverse_appear[3][k-1] + 1;
            l = count[8]-1 + reverse_appear[3][l];
        }

        if (k > l)
        {
            k = 1;
            l = N-1;
            z = z + 1;
        }

        D[i] = z;        
    }
}

SA_interval bing(SA_interval a, SA_interval b)      //并运算
{
    //将两个区间并起来
    int i;
    int num = b.index;      //b有多少个区间
    for (i = 0;i <= num;i++)
    {
        a.index += 1;    //a的index加一  表示a的区间个数多一个
        a.up[a.index] = b.up[i];
        a.down[a.index] = b.down[i];    //把b的上下界复制过来
    }

    return a;   //把合并好的集合返回

}

SA_interval InexRecur(char read[], int i, int z, int k, int l)
{
    //论文中递归代码的实现

    SA_interval this_interval;  //区间集合
    char dict[4] = {'A', 'C', 'G', 'T'};   
    int j;

    if (z < D[i])
    {
        this_interval.index = -1;   //表示没有上下界  为空
        return this_interval;
    }        
    if (i < 0)
    {
        this_interval.up[this_interval.index+1] = k;
        this_interval.down[this_interval.index+1] = l;
        this_interval.index++;
        return this_interval;
    }
    this_interval.index = -1;   //初始化为空

    //并上返回的区间 也就是将返回的区间都复制到现在的区间集合中
    this_interval = bing(this_interval, InexRecur(read, i-1, z-1, k,l));
    for (j = 0;j < 4;j++)   //分ACGT四个分支
    {
        k = count[j*2 + 2]-1 + appear[j][k-1] + 1;
        l = count[j*2 + 2]-1 + appear[j][l];
        if (k <= l)
        {            
            this_interval = bing(this_interval, InexRecur(read, i, z-1, k,l));
            if ((i==0&&read[i]=='A') || (i==1&&read[i]=='C') ||
                (i==2&&read[i]=='G') || (i==3&&read[i]=='T'))
            {
                    this_interval = bing(this_interval, InexRecur(read, i-1, z, k,l));
            }
            else
            {
                this_interval = bing(this_interval, InexRecur(read, i-1, z-1, k,l));
            }            
        }
    }
    return this_interval;
}

int find_INexact_SA_interval(int num, int file_num, int z)
{
    //对于第file_num个文件所有reads中第num个read    找对应的区间
    //写入文件   返回是否找到合法区间 1为找到   -1  没找到

    int i,j;
    int R_down,R_up;    //区间上下界  
  
    char read[70];  //一个read70个字符 

    if (file_num == 1)      //得到要比对的read
    {
        memcpy(read,reads1[num],sizeof(char)*70);   //得到要比对的文件一的第num条read 
    }
    else
    {
        memcpy(read,reads2[num],sizeof(char)*70);   //得到要比对的文件二的第num条read 
    }
    
    calculate_D(read);      //确定是哪一条read后   先计算这条read的D数组

    //接下来递归计算区间
    inEXA_intervals = InexRecur(read, 69, 1000, 1, N-1);    //这里z先指定是1000
    

    // INexact_write_file(R_up, R_down, file_num, num, 1); 
    return 1;       //成功找到对应的区间  返回1   
}

int reverse_find_INexact_SA_interval(int num, int file_num, int z)
{
    //直接在函数里写入文件   返回是否找到合法区间

    int i,j;
    int R_down,R_up;

    char orig_read[70];   //原read
    char read[70];      //反向互补后的read

    if (file_num == 1)      //得到要比对的read
    {
        //得到要比对的文件一的第num条read 
        memcpy(orig_read,reads1[num],sizeof(char)*70);
    }
    else
    {
        memcpy(orig_read,reads2[num],sizeof(char)*70);   //得到要比对的文件二的第num条read 
    }
    
    //计算这条read的反向互补
    for (j = 69;j >= 0;j--) //得到要比对的read的反向互补
    {
        if(orig_read[j] == 'A')
        {
            read[70-j-1] = 'T';
        }
        else if(orig_read[j] == 'C')
        {
            read[70-j-1] = 'G';
        }
        else if(orig_read[j] == 'G')
        {
            read[70-j-1] = 'C';
        }
        else if(orig_read[j] == 'T')
        {
            read[70-j-1] = 'A';
        }
    }    


    calculate_D(read);      //计算这条read的D数组
    
    //接下来递归计算区间
    inEXA_intervals = InexRecur(read, 69, 1000, 1, N-1);    //这里z先指定是1000

    // exact_write_file(R_up, R_down, file_num, num, -1);   
    return 1;   
}

void inexact_aligment()
{
    int i,j,flag; 
    int z;  //允许错误的上限  

    //先按文件一比对  如果找不到区间的话  就求反向互补比对
    for (i = 0;i < 10;i++)
    {
        flag = find_INexact_SA_interval(i, 1, z);  //传入read序号   先正向找
        if (flag = -1)
        {
            flag = reverse_find_INexact_SA_interval(i, 1, z);   //正向找不到  找反向互补
        }
    }

    //再比文件二  先正向找  如果不正确的话再反向互补找
    for (i = 0;i < 10;i++)
    {
        flag = find_INexact_SA_interval(i, 2, z);  //传入read序号   先正向找
        if (flag = -1)
        {
            flag = reverse_find_INexact_SA_interval(i, 2, z);   //正向找不到  找反向互补
        }
    }
}


int main()
{
    get_BWT();  //得到BWT索引
    get_count_and_O();      //得到count数组和O数组
    exact_aligment();

    get_reverse_BWT();  //得到BWT索引
    get_reverse_O();      //得到count数组和O数组
    inexact_aligment();
}
