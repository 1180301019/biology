#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#define HASHSIZE 50000
//先建一个不是全部read的图

struct vnode 
{//顶点               
    char vertex[49];      //顶点的值    是一个k-1 mer 
    struct vnode* edge[2000];//指向下一个k-1 mer
    int n;  //边用到了哪个下标
    // int index ;//在图中的顶点下标
}; 
typedef struct vnode vnode;

  
struct graph
{         //图 
    vnode vexlist[500000];   //顶点列表
    int n;                //顶点个数
};
typedef struct graph *graph;


struct __HashEntry{     //hash表里存的元素  有两个域 一个是代表的k-1mer值  另一个指向下一个元素
    char* str_value;       //即k-1mer
    int index_in_graph;      //在图中的下标值
    HashEntry *next;    //指向下一个元素  键相等的
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
char reads1[10][75];   //文件一的read
char reads2[10][75];   //文件二的read


void get_reads()
{
    FILE *fp1 = NULL;   //第一个read文件
    FILE *fp2 = NULL;   //第二个read文件    

    int i = 0,j,index = 0,len;
    char c;
    char buf1[1024];  /*缓冲区*/
    char buf2[1024];  /*缓冲区*/   


    fp1 = fopen("NC_008253_1.fastq", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式
    fp2 = fopen("NC_008253_2.fastq", "r"); //第一个逗号前是文件位置。逗号之后是打开文件方式

    if(fp1 == NULL ||fp2 == NULL)
    {
        printf("打开文件失败");
        return;
    }

    //先将文件一和文件二的read读入内存  这里读10行
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


////////////////////////////////hash//////////////////////////////////
//初始化散列表
int InitHashTable()
{
	int i;
	
	h.size=HASHSIZE;    //哈希表大小
    h.bucket = (HashEntry*)malloc(HASHSIZE*sizeof(HashEntry*));     //申请初始哈希表的桶
	
	for(i = 0;i < HASHSIZE;i++)
		h.bucket[i] = NULL;     //初始化为空指针
	return 1;
}

int hash_string(char *str)       //得到对应String值对应的整数
{
    int c,i;
    long long_hash = 0;   //映射到的桶值 即哈希表中的下标值
    int hash;

    //字符串一共49位  当成整数
    for (i = 48;i >=0;i++)
    {
        long_hash = long_hash * (49-i) * 10 + (int)str[i];
    }

    if (long_hash < 0)
        long_hash *= -1;

    hash = long_hash % 50000;

    return hash;
}

void hash_add(char* k_1_mer, int index_of_new_vert)
{
    // 向哈希表中添加元素
	// 
	int hash = hash_string(k_1_mer);    //得到k-1mer的hash值
    HashEntry* temp;

    // 新分配一个HashEntry的节点
    HashEntry* p = (HashEntry*)malloc(sizeof(HashEntry));
    memcpy(p->str_value,k_1_mer,sizeof(char)*49);   //分配字符串
    p->index_in_graph = index_of_new_vert;  //分配图中的下标

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

void construct_graph()
{
    int i,j,edge_index,left_index,right_index,index_of_new_vert;  
    g = (graph)malloc(sizeof(struct graph));        //在堆上为图申请内存
    char this_kmer[50];     
    char left_k_1_mer[49];
    char right_k_1_mer[49];

    get_reads();       //将所有的read读入内存
    

    for (i = 0;i < 10;i++)  //对于每一个read
    {
        for (j = 0;j < 70-50+1;j++)     //对于每一个k-mer   这里是50mer
        {
            memcpy(this_kmer,reads1[i][j],sizeof(char)*50);     //将对应的kmer得到
            memcpy(left_k_1_mer,this_kmer[0],sizeof(char)*49);    //得到左边的k-1mer
            memcpy(right_k_1_mer,this_kmer[1],sizeof(char)*49);   //得到右边的k-1mer

            left_index = contain(left_k_1_mer);       //如果在图中不存在返回-1  否则返回在图中的下标
            right_index = contain(right_k_1_mer);

            if (left_index == -1)     //对于左边的k-1 mer  如果图中还没有这个顶点 就加入这个顶点
            {
                g->n++; //图的顶点下标加一  
                index_of_new_vert = g->n;                  
                memcpy(g->vexlist[index_of_new_vert].vertex,left_k_1_mer,sizeof(char)*49);//把左面k-1mer加入图  
                hash_add(left_k_1_mer, index_of_new_vert);  //加入哈希表
                

                if (right_index == -1)  //如果右面k-1mer图中也没有的话 加入图 并加边
                {
                    //加顶点
                    g->n++; //图的顶点下标加一  
                    index_of_new_vert = g->n;      
                    memcpy(g->vexlist[index_of_new_vert].vertex,right_k_1_mer,sizeof(char)*49);//把右面k-1mer加入图
                    hash_add(right_k_1_mer, index_of_new_vert);  //加入哈希表

                    //加边
                    g->vexlist[g->n-1].n++;  //左边k-1mer的出边下标值加一  加上即将加上的边  出边个数加一
                    edge_index =  g->vexlist[g->n-1].n;        //得到左边k-1mer的出边下标值                    
                    g->vexlist[g->n-1].edge[edge_index] = &(g->vexlist[g->n]);  //加边  左面k-1mer指向右边                    
                } 
                else    //如果右面k-1 mer图中有  直接加边即可
                {
                    //加边  找到对应的右面k-1mer   即在图中的顶点的下标
                    g->vexlist[g->n].n++;  //左边k-1mer的出边下标值加一  加上即将加上的边  出边个数加一
                    edge_index = g->vexlist[g->n].n;    //得到左边k-1mer的出边的个数
                    g->vexlist[g->n].edge[edge_index] = &(g->vexlist[right_index]);  //加边  左面k-1mer指向右边                    
                } 
            }
            else    //如果图中已存在左面的k-1 mer
            {
                if (right_index == -1)  //如果右面k-1mer图中也没有的话 加入图 并加边
                {
                    //加顶点
                    g->n++; //图的顶点个数加一 
                    index_of_new_vert = g->n;                
                    memcpy(g->vexlist[index_of_new_vert].vertex,right_k_1_mer,sizeof(char)*49);//把右面k-1mer加入图
                    hash_add(right_k_1_mer, index_of_new_vert);  //加入哈希表

                    //加边 
                    g->vexlist[left_index].n++;
                    edge_index = g->vexlist[left_index].n;        //得到左边k-1mer的出边的个数                   
                    g->vexlist[left_index].edge[edge_index] = &(g->vexlist[g->n]);  //加边  左面k-1mer指向右边                    
                } 
                else    //如果右面k-1 mer图中有  直接加边即可
                {
                    //加边  找到对应的右面k-1mer   即在图中的顶点的下标
                    g->vexlist[left_index].n++;
                    edge_index = g->vexlist[left_index].n;        //得到左边k-1mer的出边的个数                   
                    g->vexlist[left_index].edge[edge_index] = &(g->vexlist[right_index]);  //加边  左面k-1mer指向右边 
                    
                } 
            }
        }
    }
}




int main()
{
    construct_graph();
}