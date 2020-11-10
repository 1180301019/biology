#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#define N 101   //一共5009545长  去掉头部5009476长  每行70个  一共70556行  基因组一共70*70556=4938920长 + $
#define block_num 10
#define l = 10   //前面14块的长度   4938920/15 = 329261   最后一块329266长
#define last_l = 11   //最后一块329266长  加上一个$符 长329267


typedef struct suffix{
    int SA;  //该后缀的SA
    int CSA;   //该后缀的CSA
    int index;    //该后缀的索引
    int order;
    char* str;  //该后缀 具体是串开始的地址
 }suffix;//后缀

char* T=NULL;
char* Ti[10];
int first_order;
int last_order;
suffix* T_in_order;
int length_of_T;
int count[10] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};    //分别标识$ A区  C  G  T的开始位置和结束位置   前两个一定是0 0

//获得文本串
int get_T()
{
    // char T[4938920]; //基因组
    T = (char*)malloc(120);

    FILE *fp = NULL;
    int i = 0,j;
    char c;

    fp = fopen("bio_test.fna", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
    if(fp == NULL)
    {
        printf("打开文件失败");
        return -1;
    }

    //读头
    for (i = 0;i < 69;i++)
    {
        fscanf(fp,"%c",&c);
    }

    //读基因组
    i = 0;
    //while (1)
    for (j = 0;j < 102;j++)
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
    if (i == 9)
        qsort(suf, 12, sizeof(suffix), compare_string);
    else
        qsort(suf, 10, sizeof(suffix), compare_string);
}

//得到SA的逆
void get_SA_inverse(suffix* suf, int i)
{
    if (i == 9)
        qsort(suf, 12, sizeof(suffix), compare_SA);   
    else
        qsort(suf, 10, sizeof(suffix), compare_SA);
}

//得到乱序的CSA
void out_order_CSA(suffix* suf)
{
    int i;//最后一组的长度
    int first_CSA = suf[0].index;//暂存一下第一个的index
    for (i = 0;i < 11;i++)
    {
        suf[i].CSA = suf[i+1].index;//得到CSA
    }
    suf[11].CSA = first_CSA;  //最后一个的
}


void get_CSA(suffix* suf, int i)
{
    if (i == 9)
        qsort(suf, 12, sizeof(suffix), compare_index);
    else
        qsort(suf, 10, sizeof(suffix), compare_index);
}

//提取出每一个后缀
suffix* get_suf(char* T, int i)
{
    int j;
    //首先提取出每一个后缀  然后排序
    if (i == 9)
    {
        suffix* suf = (suffix*)malloc(12*sizeof(suffix));//最后一个分组要329267块   俩*是由于数组里存的是指针
        for (j = 0;j < 12;j++)
        {
            
            suf[j].SA = j + i*10; //几后缀  
            suf[j].CSA= 0;
            suf[j].index = j;
            suf[j].order = 0;
            suf[j].str = T + j;
        }
        return suf;

    }
    else
    {
        suffix* suf = (suffix*)malloc(10*sizeof(suffix));//其他分组要329267块
        for (j = 0;j < 10;j++)
        {
            suf[j].SA = j + i*10;
            suf[j].CSA= 0;
            suf[j].index = j;
            suf[j].order = 0;
            suf[j].str = T + j;
        }
        return suf;
    }
}

void get_order(suffix* suf, int i)
{
    int this_order = -1;
    int j;
    if (i == 9)  //如果加了一个字符的串，那它依赖的order就是first_order
    {
        first_order = T_in_order[0].CSA;
        if (suf->str[0] == 'A')  //第一个字符是A，只在A区找
        {
            if (count[2] == -1)    //如果就直接没有A区  返回C区的上一个数
            {
                if (count[4] == -1)     //C区也没有   G区
                {
                    if (count[6] == -1)     //G区也没有   T区
                    {
                        if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                            suf->order = length_of_T - 1;
                        else
                            suf->order = count[8] - 1; 
                    }
                    else
                        suf->order = count[6] - 1;
                                       
                }   
                else
                    suf->order = count[4] - 1;
            } 
            else
            {
                for (j = count[2]; j <= count[3]; j++)
                {
                    if (T_in_order[j].CSA <= first_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[2] - 1;                
                else    //不为空应等于lc-1
                    suf->order = this_order;
            }
            
            

            last_order = suf->order;
        }

        else if (suf->str[0] == 'C')  //第一个字符是A，只在A区找
        {
            if (count[4] == -1)    //如果就直接没有C区  返回G区的上一个数
            {
                if (count[6] == -1)     //G区也没有   T区
                {
                    if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                        suf->order = length_of_T - 1;
                    else
                        suf->order = count[8] - 1;                    
                }   
                else
                    suf->order = count[6] - 1;
            }      
            else
            {
                for (j = count[4]; j <= count[5]; j++)
                {
                    if (T_in_order[j].CSA <= first_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[4] - 1;                
                else    //不为空应等于lc-1
                    suf->order = this_order;

                
            }

            last_order = suf->order;
            
        }

        else if (suf->str[0] == 'G')  //第一个字符是A，只在A区找
        {
            if (count[6] == -1)     //G区也没有   T区
            {
                if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                    suf->order = length_of_T - 1;
                else
                    suf->order = count[8] - 1;                    
            }   
            else
            {
                for (j = count[6]; j <= count[7]; j++)
                {
                    if (T_in_order[j].CSA <= first_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[6] - 1;                
                else    //不为空应等于lc-1
                    suf->order = this_order;
            }
            
            
            
            last_order = suf->order;
        }

        else if (suf->str[0] == 'T')  //第一个字符是A，只在A区找
        {
            if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                    suf->order = length_of_T - 1;
            else
            {
                for (j = count[8]; j < length_of_T; j++)
                {
                    if (T_in_order[j].CSA <= first_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[8] - 1;
                else    //不为空应等于lc-1
                    suf->order = this_order;

                
            }
            last_order = suf->order;
        }
    }
    else    //只改了一半
    {
        if (suf->str[0] == 'A')  //第一个字符是A，只在A区找
        {
            if (count[2] == -1)    //如果就直接没有A区  返回C区的上一个数
            {
                if (count[4] == -1)     //C区也没有   G区
                {
                    if (count[6] == -1)     //G区也没有   T区
                    {
                        if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                            suf->order = length_of_T - 1;
                        else
                            suf->order = count[8] - 1; 
                    }
                    else
                        suf->order = count[6] - 1;
                                       
                }   
                else
                    suf->order = count[4] - 1;
            } 
            else
            {
                for (j = count[2]; j <= count[3]; j++)
                {
                    if (T_in_order[j].CSA <= last_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[2] - 1;
                    
                else    //不为空应等于lc-1
                    suf->order = this_order;
            }
            last_order = suf->order;
        }

        else if (suf->str[0] == 'C')  //第一个字符是A，只在A区找
        {
            if (count[4] == -1)    //如果就直接没有C区  返回G区的上一个数
            {
                if (count[6] == -1)     //G区也没有   T区
                {
                    if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                        suf->order = length_of_T - 1;
                    else
                        suf->order = count[8] - 1;                    
                }   
                else
                    suf->order = count[6] - 1;
            } 
            else
            {
                for (j = count[4]; j <= count[5]; j++)
                {
                    if (T_in_order[j].CSA <= last_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[4] - 1;                
                else    //不为空应等于lc-1
                    suf->order = this_order;

                
            }
            last_order = suf->order;
            
        }

        else if (suf->str[0] == 'G')  //第一个字符是A，只在A区找
        {
            if (count[6] == -1)     //G区也没有   T区
            {
                if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                    suf->order = length_of_T - 1;
                else
                    suf->order = count[8] - 1;                    
            }   
            else
            {
                for (j = count[6]; j <= count[7]; j++)
                {
                    if (T_in_order[j].CSA <= last_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[6] - 1;                
                else    //不为空应等于lc-1
                    suf->order = this_order;
            }
            last_order = suf->order;
        }

        else if (suf->str[0] == 'T')  //第一个字符是A，只在A区找
        {
            if (count[8] == -1)     //T区也没有  就是只有A区，就插最后一个去
                    suf->order = length_of_T - 1;
            else
            {
                for (j = count[8]; j < length_of_T; j++)
                {
                    if (T_in_order[j].CSA <= last_order && T_in_order[j].index > this_order)
                        this_order = T_in_order[j].index;
                }
                if (this_order == -1)//集合b为空    应等于
                    suf->order = count[8] - 1;
                else    //不为空应等于lc-1
                    suf->order = this_order;
            }
            last_order = suf->order;
        }
    }
}

void update_count()
{
    int i;
    int count_$ = 0, count_A = 0, count_C = 0, count_G = 0, count_T = 0;
    for (i = 0;i < length_of_T;i++)
    {
        if (T_in_order[i].str[0] == '$')
            continue;
        else if (T_in_order[i].str[0] == 'A')
            count_A++;
        else if (T_in_order[i].str[0] == 'C')
            count_C++;
        else if (T_in_order[i].str[0] == 'G')
            count_G++;
        else if (T_in_order[i].str[0] == 'T')
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

suffix* merge_two_sequence(suffix* T, suffix* suf)
{
    suffix* new_T = (suffix*)malloc((length_of_T+10)*sizeof(suffix));   //原排好串加上新块的长度
    int i;
    for (i = 0;i < length_of_T;i++)
    {
        printf("%d\n",T[i].index);
        new_T[T[i].index] = T[i];
    }
    for (i = 0;i < 10;i++)
    {
        new_T[suf[i].index] = suf[i];
    }

    //更新T的长度
    length_of_T = length_of_T + 10;

    return new_T;
}


//分块  排最后一组 
suffix* basic_step()
{
    // long N = get_file_length();  //得到文件长度  总长5009545-头部69=5009476 每行70个 一共70556行  70*70556=4938920长 
    get_T();  //得到文本的串
    suffix* suf;    

    //切成15块  Ti存的是每块的开始位置
    int i = 0, j=0;
    for (i = 0;i < 10;i++)
    {
        Ti[i] = T + i*10;
        // Ti[i]=&(T[i*10]);
        printf("%s\n", Ti[i]);
    } 

    //排最后一块
    suf = get_suf(Ti[9], 9);  //得到所有后缀
    sort_suf(suf, 9);         //得到SA
    //调整下标
    for (i = 0;i < 12;i++)
    {
        suf[i].index = i;
    }
    get_SA_inverse(suf, 9);      //得到SA的逆
    first_order = suf[0].index;
    out_order_CSA(suf);
    get_CSA(suf, 9);

    length_of_T = 12;    ///更新T的长度  这里应该等于最后一组的长度

    return suf;
}

suffix* merge_step()
{
    suffix* suf; 
    int i, j, block_i;
    for (block_i = 8; block_i >= 0;block_i--)
    {

    
    ////////////////////////////////////////////////step1  排序///////////////////////////////////////////////////////
    suf = get_suf(Ti[block_i], block_i);
    sort_suf(suf, block_i);         //     排序
    //调整下标
    for (i = 0;i < 10;i++)
    {
        suf[i].index = i;
    }
    get_SA_inverse(suf, block_i);      //  让这个块的顺序是依次递减一个字符的

    //////////////////////////////////////////////step2     得到order///////////////////////////////////////////////////////
    for(i = 9;i >=0;i--)        //对块内的每一个后缀获得order
    {
        get_order(&suf[i], i);
    }

    //////////////////////////////////////////////step3     融合///////////////////////////////////////////////////////

    int* f = (int *)malloc(length_of_T*sizeof(int));
    int* g = (int *)malloc(10*sizeof(int));    //一组的大小
    
    
    //重新计算index   分别计算原串的和新来块的
    //原串的
    int T_index = 0, suf_index = 0;
    for (i = 0;i < length_of_T;i++)
    {
        for (j = 0;j < 10;j++)  //统计suf中排在该串前面的
        {
            if (suf[j].order < T_in_order[i].index)
                T_index++;
        }
        T_in_order[i].index = T_index + T_in_order[i].index;    //原本的排名加上插在它前面的
        f[i] = T_in_order[i].index;     //记在f数组中
        T_index = 0;       //每次计算后重置以下
    }

    //新块的
    for (j = 0;j < 10;j++)
    {
        suf[j].index = suf[j].index + suf[j].order + 1;     //新块
        g[j] = suf[j].index;          //在g数组中记一下
    }



    //合并为一个块
    T_in_order = merge_two_sequence(T_in_order, suf);

    // 更新CSA值
    int jf = 1, jg = 0, t, temp;
    suffix ttt;
    int CSA_0 = f[T_in_order[0].CSA];   //提前记下排好串的最长的新位置   便于给后面新块最后一个的CSA赋值
    // for (t = 1;t < length_of_T-1;t++)
    // {
    //     if (jg == 9)   //新块最短那个  就相当于刚加了一个字符的串   它的CSA对应得串在上一块  单独赋值
    //     {
    //         T_in_order[g[jg]].CSA = CSA_0;//最长串所在位置
    //     }
    //     if (t == f[jf])     //如果是原串的话  就是去找它对应的CSA那个串的新位置
    //     {
    //         ttt = T_in_order[t];
    //         T_in_order[t].CSA = f[T_in_order[t].CSA];
    //         jf++;
    //     }
    //     else    //如果是新块的话    就是依次找  去头
    //     {
    //         if (jg == 0)        //最长串的位置赋给$的CSA
    //         {
    //             T_in_order[0].CSA = g[jg];
    //         }
    //         suffix temm = T_in_order[t];
    //         T_in_order[g[jg]].CSA = g[jg+1];
    //         jg++;
    //     }        
    // }    
            
    for (t = 1;t < length_of_T;t++)   //length_of_T-2次    f中即原串减一 是减掉了$符   G中减一 是减掉了最后的
    {
        
        if (t == f[jf])     //如果是原串的话  就是去找它对应的CSA那个串的新位置
        {
            ttt = T_in_order[t];
            T_in_order[t].CSA = f[T_in_order[t].CSA];
            jf++;
        }
        else    //如果是新块的话    就是依次找  去头
        {
            
            suffix temm = T_in_order[t];
            T_in_order[g[jg]].CSA = g[jg+1];
            jg++;
        }        
    }
    //新块最短那个  就相当于刚加了一个字符的串   它的CSA对应得串在上一块  单独赋值
    T_in_order[g[9]].CSA = CSA_0;//最长串所在位置
    //最长串的位置赋给$的CSA
    T_in_order[0].CSA = g[0];

    //更新count
    update_count();
    }



    printf("0");
}

int main()
{
    T_in_order = basic_step();
    update_count();

    suffix* final_suf = merge_step();


    
    return 0;
}

