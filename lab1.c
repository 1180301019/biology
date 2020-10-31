#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#define N 4938921   //一共5009545长  去掉头部5009476长  每行70个  一共70556行  基因组一共70*70556=4938920长 + $
#define block_num 15
#define l = 329261   //前面14块的长度   4938920/15 = 329261   最后一块329266长
#define last_l = 329267   //最后一块329266长  加上一个$符 长329267


char* T=NULL;
char* Ti[15];



typedef struct suffix{
    int SA;  //该后缀的SA
    int CSA;   //该后缀的CSA
    int index;    //该后缀的索引
    char* str;  //该后缀 具体是串开始的地址
 }suffix;//后缀




// long get_file_length()
// {
//     FILE *fp = NULL;

//     fp = fopen("NC_008253.fna", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
//     if(fp == NULL)
//     {
//         printf("打开文件失败");
//         return -1;
//     }
//     fseek(fp,0,SEEK_END);//定位到文件的最后面
//     long text_length  = ftell(fp);//ftell获得该文件指示符此时的偏移量,此时已经是在文件末尾,故能获得文件的大小
//     printf("%ld\n",text_length);

//     fp = fopen("NChead", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
//     if(fp == NULL)
//     {
//         printf("打开文件失败");
//         return -1;
//     }
//     fseek(fp,0,SEEK_END);//定位到文件的最后面
//     long head_length  = ftell(fp);//ftell获得该文件指示符此时的偏移量,此时已经是在文件末尾,故能获得文件的大小
//     printf("%ld\n",head_length);

//     long length = text_length - head_length;

//     fclose(fp);  //记得用完关闭文件
//     return length;
// }

//获得文本串
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
    //while (1)
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

// /////////////////////////////////////////////////获得SA/////////////////////////////////////////////////////////

//比较两个串的大小
int compare_string(const void * a, const void * b)  
{
    //return strcmp(((suffix*)a)->str, ((suffix*)b)->str);    
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

// ///////////////////////////////////////////////获得SA的逆 就是得到后缀排的index////////////////////////////////////////////////////
int compare_SA(const void * a,const void * b)
{
    return ((suffix*)a)->SA > ((suffix*)b)->SA? 1:-1;
}

//////////////////////////////////////////得到正序的CSA   就是排index////////////////////////////////////////////////
//比较两个串的大小
int compare_index(const void * a,const void * b)
{
    return ((suffix*)a)->index > ((suffix*)b)->index? 1:-1;
}

/////////////////////////////////////////////////////////////////////////////////////////////
void sort_suf(suffix* suf, int i)
{
    if (i == 14)
        qsort(suf, 329267, sizeof(suffix), compare_string);
    else
        qsort(suf, 329261, sizeof(suffix), compare_string);
}

//得到SA的逆
void get_SA_inverse(suffix* suf)
{
    qsort(suf, 329267, sizeof(suffix), compare_SA);   
}

//得到乱序的CSA
void out_order_CSA(suffix* suf)
{
    int i;//最后一组的长度
    int first_CSA = suf[0].index;//暂存一下第一个的index
    for (i = 0;i < 329266;i++)
    {
        suf[i].CSA = suf[i+1].index;//得到CSA
    }
    suf[329266].CSA = first_CSA;  //最后一个的
}

void get_CSA(suffix* suf)
{
    qsort(suf, 329267, sizeof(suffix), compare_index);
}

//提取出每一个后缀
suffix* get_suf(char* T, int i)
{
    int j;
    //首先提取出每一个后缀  然后排序
    if (i == 14)
    {
        suffix* suf = (suffix*)malloc(329267*sizeof(suffix));//最后一个分组要329267块   俩*是由于数组里存的是指针
        for (j = 0;j < 329267;j++)
        {
            
            suf[j].SA = j + i*14; //几后缀  
            suf[j].CSA= 0;
            suf[j].index = j;
            suf[j].str = T + j;
        }
        return suf;

    }
    else
    {
        suffix* suf = (suffix*)malloc(329261*sizeof(suffix));//其他分组要329267块
        for (j = 0;j < 329261;j++)
        {
            suf[j].SA = j + i*10;
            suf[j].CSA= 0;
            suf[j].index = j;
            suf[j].str = T + j;
        }
        return suf;
    }
}



//分块  排最后一组 
suffix* basic_step()
{
    // long N = get_file_length();  //得到文件长度  总长5009545-头部69=5009476 每行70个 一共70556行  70*70556=4938920长 
    get_T();  //得到文本的串
    suffix* suf14;

    //切成15块  Ti存的是每块的开始位置
    char* Ti[15] ;
    int i = 0, j=0;
    for (i = 0;i < 15;i++)
    {
        // Ti[i] = T + i*329261;
        Ti[i]=&(T[i*329261]);
    } 



    //排最后一块
    suf14 = get_suf(Ti[14], 14);  //得到所有后缀
    sort_suf(suf14, 14);         //得到SA
    //调整下标
    for (i = 0;i < 329267;i++)
    {
        suf14[i].index = i;
    }
    get_SA_inverse(suf14);      //得到SA的逆
    out_order_CSA(suf14);
    get_CSA(suf14);

    return suf14;


    // 验证一下地址是不是正确
    // for (i = 0;i < 15;i++)
    // {
    //     printf("%d\n",strlen(Ti[i]));
    // }

    
    





    // printf("%s\n",T);











//    printf("%ld\n",block_num);
//    printf("%ld\n",l);




    //test
//    FILE *fp = NULL;
//    int i = 0;
//    char c;
//
//    fp = fopen("NC_008253.fna", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
//    if(fp == NULL)
//    {
//        printf("打开文件失败");
//        return -1;
//    }
//    for (i = 0;i < 70;i++)
//    {
//        if (fscanf(fp,"%c",&c)!=EOF)
//            printf("%c",c);
//    }
//
//    fclose(fp);

    printf("0");


}

int main()
{
    basic_step();
    return 0;
}

