#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#define HASHSIZE 50000
#define GRAPHSIZE 64000000
//每个read150个字符       640000个read

struct vnode 
{//顶点               
    char vertex[49];      //顶点的值    是一个k-1 mer 
    int edge[2000];
    // struct vnode* edge[2000];//指向下一个k-1 mer
    int n;  //边用到了哪个下标
    int index ;//在图中的顶点下标
}; 
typedef struct vnode vnode;

  
struct graph
{         //图 
    int n;                //顶点个数
    int size;
    vnode **vexlist;       //顶点链表
};
typedef struct graph *graph;


struct __HashEntry{     //hash表里存的元素  有两个域 一个是代表的k-1mer值  另一个指向下一个元素
    char* str_value;       //即k-1mer
    int index_in_graph;      //在图中的下标值
    struct __HashEntry *next;    //指向下一个元素  键相等的
};
typedef struct __HashEntry HashEntry;//哈希表所保存元素（即键值对 《key, value》）类型

struct __HashTable{
    HashEntry **bucket;        
    int size;           //哈希表的大小
};
typedef struct __HashTable HashTable; //哈希表，其中 bucket 指向大小为size的、元素类型为 HashEntry*的指针数组
                                        //链地址法处理冲突
//全局变量部分
graph g;    //图
HashTable h;
char reads1[640000][150];   //文件一的read
// char reads2[640000][150];   //文件二的read


void get_reads()
{
    FILE *fp1 = NULL;   //第一个read文件   

    int i = 0,j,index = 0,len;
    char c;
    char buf1[1024];  /*缓冲区*/  


    fp1 = fopen("NC_008253_1.fastq", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式

    if(fp1 == NULL)
    {
        printf("打开文件失败");
        return;
    }

    //先将文件一和文件二的read读入内存  这里除以4是read数
    for (i = 0;i < 640000*4;i++)
    {        
        fgets(buf1,1024,fp1);      //按行读取文件1
        len = strlen(buf1);         //得到本行的长度
        buf1[len-1] = '\0';  //去掉换行符
        if (i % 4 != 1)     //观察文件可知，第二行是reads，每4行是一个reads
            continue;
        else
        {
            memcpy(reads1[index],buf1,sizeof(char)*150);     //将缓冲区的read写入数组
            index++;    //索引加一
        }
    }

    fclose(fp1);
}


////////////////////////////////hash//////////////////////////////////
//初始化散列表
int InitHashTable()
{
	int i;
	
	h.size=HASHSIZE;    //哈希表大小
    h.bucket = (HashEntry**)malloc(HASHSIZE*sizeof(HashEntry*));     //申请初始哈希表的桶
	
	for(i = 0;i < HASHSIZE;i++)
		h.bucket[i] = NULL;     //初始化为空指针
	return 1;
}

int hash_string(char *str)       //得到对应String值对应的整数
{
    int c,i;
    long long_hash = 0;   //映射到的桶值 即哈希表中的下标值
    int hash;

    for (i = 48;i >=0;i--)
    {
        long_hash = (long_hash*4 + (int)str[i]) % 50000;
    }   
    hash = long_hash % 50000;

    return hash;
}

void hash_add(char* k_1_mer, int index_of_new_vert)
{
    // 向哈希表中添加元素
	
	int hash = hash_string(k_1_mer);    //得到k-1mer的hash值
    HashEntry* temp;

    // 新分配一个HashEntry的节点
    HashEntry* p = (HashEntry*)malloc(sizeof(HashEntry));
    // memcpy(p->str_value,k_1_mer,sizeof(char)*50);   //分配字符串
    p->str_value = k_1_mer;
    p->index_in_graph = index_of_new_vert;  //分配图中的下标
    p->next = NULL;

    if (h.bucket[hash] == NULL)     //如果对应的桶内还没有值的话
    {
        h.bucket[hash] = p;
    }
    else    //桶内有值 链到前面即可
    {
        temp = h.bucket[hash];
        h.bucket[hash] = p;
        p->next = temp;
    }
}

int contain(char* k_1_mer)
{
    int hash = hash_string(k_1_mer);    //得到k-1mer的hash值
    HashEntry* temp;
    int cmp_result;

    if (h.bucket[hash] == NULL)
    {
        return -1;  // 该k-1mer在hash表中不存在
    }
    else
    {
        temp = h.bucket[hash];
        while (temp != NULL)
        {
            cmp_result = strcmp(k_1_mer, temp->str_value);
            if (cmp_result == 0)
            {
                return temp->index_in_graph;    //该k-1mer在图中的下标
            }    
            else
            {
                temp = temp->next;     //把桶内所有元素比完
            }            
        }
    }

    return -1;  // 该k-1mer在hash表中不存在
}
/////////////////////////////////////////////////////////////////////////

void initgraph()
{
    int i;

    g = (graph)malloc(sizeof(struct graph));        //在堆上为图申请内存

    g->size = GRAPHSIZE;
    g->vexlist = (vnode**)malloc((g->size)*sizeof(vnode*));
    g->n = -1;

    for (i = 0;i < GRAPHSIZE;i++)
        g->vexlist[i] = NULL;
    return;
}

void add_ver(char* k_1_mer)
{
    int index_of_new_vert;
    vnode *v;

    g->n++; //图的顶点下标加一  
    index_of_new_vert = g->n; 

    //先申请一个节点
    v = (vnode *)malloc(sizeof(vnode));    
    memcpy(v->vertex,k_1_mer,sizeof(char)*50);//把左面k-1mer加入图  
    v->n = -1;
    v->index = index_of_new_vert;
    
    g->vexlist[index_of_new_vert] = v;
    hash_add(k_1_mer, index_of_new_vert);  //加入哈希表
}

void construct_graph()
{
    int i,j,edge_index,left_index,right_index,index_of_new_vert;  
    char* this_kmer;     
    char* left_k_1_mer;
    char* right_k_1_mer;

    get_reads();       //将所有的read读入内存

    InitHashTable();    //初始化哈希表
    initgraph();    //初始化图
    

    for (i = 0;i < 10;i++)  //对于每一个read
    {
        for (j = 0;j < 150-50+1;j++)     //对于每一个k-mer   这里是50mer
        {
            this_kmer = &(reads1[i][j]);    //将对应的kmer得到
            left_k_1_mer = &(this_kmer[0]);     //得到左边的k-1mer
            right_k_1_mer = &(this_kmer[1]);    //得到右边的k-1mer  

            left_index = contain(left_k_1_mer);       //如果在图中不存在返回-1  否则返回在图中的下标
            right_index = contain(right_k_1_mer);

            if (left_index == -1)     //对于左边的k-1 mer  如果图中还没有这个顶点 就加入这个顶点
            {
                //把左面k-1mer加入图 
                add_ver(left_k_1_mer);                

                if (right_index == -1)  //如果右面k-1 mer图中也没有的话 加入图 并加边
                {
                    //加顶点    //把右面k-1mer加入图
                    add_ver(right_k_1_mer);

                    //加边
                    g->vexlist[g->n-1]->n++;  //左边k-1mer的出边下标值加一  加上即将加上的边  出边个数加一
                    edge_index =  g->vexlist[g->n-1]->n;        //得到左边k-1mer的出边下标值                     
                    g->vexlist[g->n-1]->edge[edge_index] = g->n;                    
                } 
                else    //如果右面k-1 mer图中有  直接加边即可
                {
                    //加边  找到对应的右面k-1mer   即在图中的顶点的下标
                    g->vexlist[g->n]->n++;  //左边k-1mer的出边下标值加一  加上即将加上的边  出边个数加一
                    edge_index = g->vexlist[g->n]->n;    //得到左边k-1mer的出边的个数
                    // g->vexlist[g->n].edge[edge_index] = &(g->vexlist[right_index]);  //加边  左面k-1mer指向右边
                    g->vexlist[g->n]->edge[edge_index] = right_index;                    
                } 
            }
            else    //如果图中已存在左面的k-1 mer
            {
                if (right_index == -1)  //如果右面k-1mer图中也没有的话 加入图 并加边
                {
                    //加顶点    //把右面k-1mer加入图
                    add_ver(right_k_1_mer);

                    //加边 
                    g->vexlist[left_index]->n++;
                    edge_index = g->vexlist[left_index]->n;        //得到左边k-1mer的出边的个数                   
                    // g->vexlist[left_index].edge[edge_index] = &(g->vexlist[g->n]);  //加边  左面k-1mer指向右边 
                    g->vexlist[left_index]->edge[edge_index] = g->n;                   
                } 
                else    //如果右面k-1 mer图中有  直接加边即可
                {
                    //加边  找到对应的右面k-1mer   即在图中的顶点的下标
                    g->vexlist[left_index]->n++;
                    edge_index = g->vexlist[left_index]->n;        //得到左边k-1mer的出边的个数                   
                    // g->vexlist[left_index].edge[edge_index] = &(g->vexlist[right_index]);  //加边  左面k-1mer指向右边 
                    g->vexlist[left_index]->edge[edge_index] = right_index;
                } 
            }
        }
    }
}

void find_contig()
{
    int i = 0,index = -1,flag = 1;  //flag是标志是否为开始串
    int index_not_one[500] = {0};
    FILE *fp1 = NULL;   //第一个read文件
    

    fp1 = fopen("contig.txt", "w"); //第一个逗号前是文件位置。逗号之后是打开文件方式

    if(fp1 == NULL)
    {
        printf("打开文件失败");
        return;
    }

    while (i < g->n)    //按顺序遍历图上的所有顶点
    {
        if (i == 0) //如果是图中第一个顶点  是开始串
        {
            if (g->vexlist[i]->n == 0)   //指的是边的下标用到0 即有一个出边
            {
                fputs(g->vexlist[i]->vertex, fp1);
                printf("%s", g->vexlist[i]->vertex);     //是开始串所以要打印一整个串
                //判断出边
                if (g->vexlist[i]->edge[0] != i+1)   //如果只有一个出边而且出边还不是按顺序的下一个顶点 那么说明这是一个k-mer的结束
                {
                    index++;
                    index_not_one[index] = i+1;     //记住这个顶点 遍历图结束后重新找这个节点
                    flag = 0;   ////表示不是开始串
                    i = g->vexlist[i]->edge[0];  //去下一个顶点
                }
                else    //有一个出边而且出边是下一个节点
                {
                    flag = 0;//表示不是开始串
                    i++;    //下标加一 去下一个顶点
                }                
            }
            else
            {
                i++;    //开始的顶点不只一个出边  忽略 去下一个顶点  并且flag=1不变 表示是开始串 需打印全串 
            }                
        }
        else    //不是图中第一个顶点
        {
            if (g->vexlist[i]->n == 0)   //指的是边的下标用到0 即有一个出边
            {
                if (flag == 1)
                {
                    fprintf(fp1, "\n");
                    fputs(g->vexlist[i]->vertex, fp1);
                    printf("\n%s", g->vexlist[i]->vertex);     //是开始串所以要打印一整个串             
                }
                else
                {
                    fputc(g->vexlist[i]->vertex[48], fp1);
                    printf("%c", g->vexlist[i]->vertex[48]);     //不是开始串  所以打印最后一个字符即可
                }

                //判断出边
                if (g->vexlist[i]->edge[0] != i+1)   //如果只有一个出边而且出边还不是按顺序的下一个顶点 那么说明这是一个k-mer的结束
                {
                    index++;
                    index_not_one[index] = i+1;     //记住这个顶点 遍历图结束后重新找这个节点
                    flag = 0;   ////表示不是开始串
                    i = g->vexlist[i]->edge[0];  //去下一个顶点
                }
                else    //有一个出边而且出边是下一个节点
                {
                    flag = 0;//表示不是开始串
                    i++;    //下标加一 去下一个顶点
                }  
            }
            else    //出边不是1 忽略该点
            {
                fputc(g->vexlist[i]->vertex[48], fp1);
                printf("%c", g->vexlist[i]->vertex[48]);
                flag = 1;//表示不是开始串
                i++;    //下标加一 去下一个顶点
            }
            
        }
    }

    fclose(fp1);
    return;
}


int main()
{
    construct_graph();
    find_contig();

    printf("0");
}