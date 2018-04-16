#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cassert>
//#include "random.h"
#include "recombination.h"
using namespace std; //应用于全局，不用管他直接用
double lowBound = 0;//自变量取值范围下限
double uppBound = 1;//自变量取值范围上限
int nvar = 30;//自变量x的个数
int nobj = 2;//目标函数y的个数
int pops = 100;//种群个数
int max_gen = 500;//运行的最大代数，可以不设成500，因为主函数有赋值
int seed = 109;//随机种子数，此值可以随意设
long rnd_uni_init;//64位的有符号整形，这个数是为了搭配random.h中的rnd_uni函数生成0到1之间的随机数
char filename[1024];
FILE *fp1, *fp2;
//实数交叉SBX里面的参数
int etax = 20;
int etam = 20;
double realx, realm; // SBX交叉概率 and 多项式变异
double realb = 0.9; // 从种群中选择供交配的0
int gID;
void execute(char *algorithm);
void minfastsort(double x[], int idx[], int n, int m);
void objective(vector <double> &x_var, vector <double> &y_obj);
double rnd_uni(long *l);
struct point
{
	float x1, x2;
	float y1, y2;
};
bool cmp(point p1, point p2)
{
	if (p1.x1 < p2.x1) return true;
	return false;
}
class zc_NSGA2Ind
{
public:
	vector <double> x_var;//动态数组，用于存放自变量x
	vector <double> y_obj;//动态数组，用于存放目标y
	int np;//被支配的数量
	int sp;//支配别人的数量
	void rnd_init();//随机初始化
	void obj_eval();//计算目标函数
	void show_objective();//此函数是显示目标函数，但是在调用过程中没用到，调试的时候可以用来检测
	zc_NSGA2Ind();//构造函数，调用此类时自动调用此函数，用于初始化
	~zc_NSGA2Ind();//析构函数，最后调用此函数，起到清零的作用
	bool operator < (const zc_NSGA2Ind &ind2);//对"<"重载,重新赋予其意义，使"<"能在个体之间进行比较
	bool operator == (const zc_NSGA2Ind &ind2);//同上，对"=="重载，
	void operator = (const zc_NSGA2Ind &ind2);//同上
};

//构造函数，调用此类时自动调用此函数，用于初始化
zc_NSGA2Ind::zc_NSGA2Ind()
{
	for (int i = 0; i < nvar; i++)
	{
		x_var.push_back(0.0);
	}
	for (int j = 0; j < nobj; j++)
	{
		y_obj.push_back(0.0);
	}
}

//析构函数，最后调用此函数，起到清零的作用
zc_NSGA2Ind::~zc_NSGA2Ind()
{
	x_var.clear();
	y_obj.clear();
}

class zc_NSGA2
{
public:
	zc_NSGA2();//构造函数，在本程序中没什么大作用
	virtual ~zc_NSGA2();//析构函数，本程序中没什么大作用
	void execute(int run);//执行函数，主函数通过此函数调用本类中的各个函数
	void init_population();//初始化个体，包括对自变量和目标函数的初始化
	int tour_selection();//轮盘赌选择个体
	void fill_union(zc_NSGA2Ind &ind);//将N个个体放到offspring中，并保证没有重复
	void rank_population();//非支配排序
	void evaluation_density();//拥挤度算子选择，此过程在非支配排序之后
	void evolution();//进化生成子代的过程，包括交叉和变异
	void save_front(char savefilename[1024]);//将最后的目标函数值输出到POF文件夹中
	void save_ps(char savefilename[1024]);//将最后的目标函数所对应的自变量输出到POS文件夹中
	vector <zc_NSGA2Ind> population;//种群个体P集合，注意此数组存放的个体类型是CNSGA2Ind
	vector <zc_NSGA2Ind> offspring;//种群P集合进化后的子代Q集合，注意同上
};

//构造函数
zc_NSGA2::zc_NSGA2() {}

//析构函数
zc_NSGA2::~zc_NSGA2() {}

//对"<"重载,重新赋予其意义，使"<"能在个体之间进行比较
bool zc_NSGA2Ind::operator < (const zc_NSGA2Ind &ind2)
{
	int flag2 = 0;
	for (int n = 0; n < nobj; n++)
	{
		if (ind2.y_obj[n] < y_obj[n])
			return false;
		if (ind2.y_obj[n] == y_obj[n])
			flag2++;
	}
	if (flag2 == nobj)
		return false;
	return true;
}

//对"=="重载，重新赋予其意思，使"=="能在个体之间进行比较
bool zc_NSGA2Ind::operator ==(const zc_NSGA2Ind &ind2)
{
	int flag = 0;
	for (int n = 0; n < nobj; n++)//nobj是目标函数的个数
	{
		if (ind2.y_obj[n] != y_obj[n])
			return false;
	}
	return true;
}

//对"="重载，重新赋予其意思，使"="能用在个体之间复制
void zc_NSGA2Ind::operator = (const zc_NSGA2Ind &ind2)
{
	for (int n = 0; n <nvar; n++)
		x_var[n] = ind2.x_var[n];
	for (int n = 0; n < nobj; n++)
		y_obj[n] = ind2.y_obj[n];
	np = ind2.np;
}

int main()
{
	srand(time(NULL));
	execute("NSGA-II");
	return 0;
}

void execute(char *algorithm)
{
	int run = 1;
	zc_NSGA2 NSGA2;
	NSGA2.execute(run);
}

void zc_NSGA2::execute(int run)
{
	seed = (seed + 23) % 1377;//设置一个随机种子
	rnd_uni_init = -(long)seed;//由种子产生一个随机数
	FILE *fp1, *fp2; /*定义文件指针*/
	FILE *fp;
	int gen = 1;//第一代
	cout << "***********第" << gen << "代!**********" << endl;
	init_population();//调用初始化函数
	cout << endl;
	for (gen = 2; gen <= max_gen; gen++)
	{
		cout << "***********第" << gen << "代!**********" << endl;
		evolution();//进化操作：交叉，变异
		rank_population();//排序
		evaluation_density();//拥挤度比较
		if (gen % max_gen == 0)
		{
			cout << "朱晨1";
			fopen_s(&fp1, "test1.txt", "w+"); /*建立一个文字文件只写*/
			fopen_s(&fp2, "test2.txt", "w+"); /*建立一个文字文件只写*/
			fopen_s(&fp, "data.txt", "w+");
			cout << "朱晨2";
			for (int n = 0; n < pops; n++)
			{
				for (int k = 0; k < nobj; k++)
				{
					fprintf(fp1, "%f ", population[n].y_obj[k]); /*向所建文件写一字符串*/
					cout << "POF = " << population[n].y_obj[k] << " ";
					cout << "朱晨3 : ";
				}
				fprintf(fp1, "\n");
			}
			for (int m = 0; m < pops; m++)
			{
				for (int s = 0; s < nvar; s++)
				{
					fprintf(fp2, "%f ", population[m].x_var[s]); /*向所建文件写一字符串*/
					cout << "POS = " << population[m].x_var[s] << " ";
					cout << "朱晨4 : ";
				}
				fprintf(fp2, "\n");
			}
			for (int n = 0; n < pops; ++n)
			{
				//for (int s = 0; s < nvar; ++s)
				//{
					//fprintf(fp, "%f ", population[n].x_var[s]);
				//}
				for (int s = 0; s < nobj; ++s)
				{
					fprintf(fp, "%f ", population[n].y_obj[s]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp1); fclose(fp2);
			fclose(fp);

			/*ifstream in("data.txt");
			vector<point> ps;
			point p;
			while (in >> p.x1 >> p.x2 >> p.y1 >> p.y2)
			{
				ps.push_back(p);
			}
			sort(ps.begin(), ps.end(), cmp);
			in.close();
			ofstream out("data.txt");
			for (int i = 0; i < ps.size(); ++i)
			{
				out << ps[i].x1 << " " << ps[i].x2 << " " << ps[i].y1 << " " << ps[i].y2 << "\n";
			}
			out.close();*/
		}
	}
	population.clear();
	offspring.clear();
}

//初始化个体，包括对自变量和目标函数的初始化
void zc_NSGA2::init_population()
{
	for (int n = 0; n < pops; n++)
	{
		zc_NSGA2Ind ind;//定义了一个类变量
		ind.rnd_init();
		population.push_back(ind);
	}
}

//随机初始化
void zc_NSGA2Ind::rnd_init()
{
	for (int n = 0; n < nvar; n++)//nvar为自变量个数，此处nvar=2
	{
		x_var[n] = lowBound + rnd_uni(&rnd_uni_init) * (uppBound - lowBound);
		printf("第%d个变量是：", n);
		std::cout << x_var[n] << " ";
		std::cout << "\n";
	}
	obj_eval();
	show_objective();
}

//计算目标函数
void zc_NSGA2Ind::obj_eval()
{
	objective(x_var, y_obj);
}

//目标函数
void objective(vector <double> &x_var, vector <double> &y_obj)
{
	//y_obj[0] = (x_var[0] * x_var[0]) - x_var[1];
	//y_obj[1] = (x_var[0] + x_var[1] - 2)*x_var[0];

	y_obj[0] = x_var[0];
	double g = 0;
	for (int i = 1; i < nvar; ++i)
	{
		g = g + x_var[i];
	}
	g = g / (nvar - 1);
	g = 1 + 9 * g;
	double h = 1 - sqrt(x_var[0] / g);
	y_obj[1] = g * h;
}

//此函数是显示目标函数
void zc_NSGA2Ind::show_objective()
{
	for (int n = 0; n < nobj; n++)//nobj是目标函数的个数
	{
		std::cout << "目标: " << n << " ";
		std::cout << y_obj[n] << "\n";
	}
}

//进化生成子代的过程，包括交叉和变异
void zc_NSGA2::evolution()
{
	for (int n = 0; n < pops; n++)
	{
		int p1, p2;
		p1 = tour_selection();//将轮盘赌选择的返回值赋给P1
		while (1)//此循环的目的在于选出两个不同的个体，以便后面进行交叉
		{
			p2 = tour_selection();//将轮盘赌选择的返回值赋给P2
			if (p1 != p2)
				break;
		}
		zc_NSGA2Ind child;
		//调用交叉函数(库函数)
		real_sbx_xover2(population[p1], population[p2], child);
		//调用变异函数(库函数)
		realmutation(child, 1.0 / nvar);
		//计算子代目标值
		child.obj_eval();
		fill_union(child);
		fill_union(population[n]);
	}
}

int zc_NSGA2::tour_selection()
{
	int p1 = int(rnd_uni(&rnd_uni_init) * pops) % pops;//产生一个随机数
	int p2 = int(rnd_uni(&rnd_uni_init) * pops) % pops;//产生第二个随机数
	if (population[p1].np < population[p2].np)
		return p1;
	else
		return p2;
}

//非支配排序，即给offspring集合中的每一个个体一个rank值
void zc_NSGA2::rank_population()
{
	int size = offspring.size();//size = 2N
	int** cset;//定义一个二维数组
	cset = new int*[size];
	for (int i = 0; i < size; i++)
		cset[i] = new int[size];
	//cset[2N][2N],cset[0][7]=96表示第0个点支配的第7个点是96
	int* rank = new int[size];
	//rank[2N]，表示被支配的点，如rank[0]=2表示第0个点被两个点支配
	cout << "size= " << size << "\n\n";
	for (int i = 0; i < size; i++)
	{
		rank[i] = 0;//支配第i个点的个数为0
		offspring[i].np = -1;//支配第i个点的个数为-1
		offspring[i].sp = 0;//第i个点所支配的集合为空
	}//初始化
	for (int k = 0; k < size; k++)                           //这里貌似可以改进一下   xct
	{
		for (int j = 0; j < size; j++)
		{
			if (k != j)
			{
				if (offspring[j] < offspring[k] && !(offspring[j] == offspring[k]))//j支配k
					rank[k]++;//支配第k个点的个数+1
				if (offspring[k] < offspring[j] && !(offspring[k] == offspring[j]))//k支配j
				{
					offspring[k].sp++;//k所支配的集合+1
					int m = offspring[k].sp - 1;
					cset[k][m] = j;//cset[2N][2N],cset[k][m]=j表示第k个点所支配的第m个点是j
				}
			}
		}
		printf("支配第%d个点的个数为%d：\n", k, rank[k]);
		printf("点%d所支配的点为：", k);
		cout << offspring[k].sp << "\n";
	}
	int curr_rank = 0;
	while (1)
	{
		int stop_count = 0;
		int* rank2 = new int[size];
		for (int k = 0; k < size; k++)
			rank2[k] = rank[k];//复制
		for (int k = 0; k < size; k++)
		{
			if (offspring[k].np == -1 && rank[k] == 0)
			{
				offspring[k].np = curr_rank; //被支配点数为0的点的排名，0为最好
				for (int j = 0; j < offspring[k].sp; j++)//在k所支配的集合内
				{
					int id = cset[k][j];//将点k所支配的第j+1个点给id
					rank2[id]--; //将支配id的点得个数减一
					stop_count++;
				}
			}
		}
		for (int k = 0; k < size; k++)
		{
			rank[k] = rank2[k];//将更新后支配k的点数复制回rank[]
		}
		delete[] rank2;
		curr_rank++;
		//cout<<"curr_rank: "<<curr_rank<<"\n\n\n";
		if (stop_count == 0)
			break;
	}
	delete[] rank;
	for (int i = 0; i < size; i++)
		delete cset[i];
	delete[] cset;
}

//拥挤度算子选择，此过程在非支配排序之后
void zc_NSGA2::evaluation_density()
{
	population.clear();
	int size = offspring.size();     //size = 2N
	int rank = 0;
	while (1)
	{
		int count = 0;
		for (int i = 0; i < size; i++)
		{
			if (offspring[i].np == rank)
				count++;
		}
		//cout<<"count: "<<count<<"\n\n";
		int size2 = population.size() + count; //最优点数量
		//cout<<"最优点数量size2: "<<size2<<"\n\n";
		if (size2 > pops)
			break;
		for (int i = 0; i < size; i++)
			if (offspring[i].np == rank)
				population.push_back(offspring[i]);//遍历所有个体，找到不受支配的点
		rank++;
		//cout<<"rank:"<<rank<<"\n\n";
		if (population.size() >= pops)
			break;
	}//将非支配解按优秀程度不断填充到population中，直到大小等于pops容量为止,如果此时还不够则按照拥挤度距离继续填充
	if (population.size() < pops)
	{
		vector <zc_NSGA2Ind> list;//保存可能溢出的那一级的个体
		for (int i = 0; i < size; i++)
			if (offspring[i].np == rank)
				list.push_back(offspring[i]);
		int s2 = list.size(); //将超过pops的那个排名的点的数量
		//cout<<"超过pops的那个排名的点的数量s2: "<<s2<<"\n\n";
		double *density = new double[s2];
		int *idx = new int[s2];
		int *idd = new int[s2];
		double *obj = new double[s2];
		for (int i = 0; i < s2; i++)
		{
			density[i] = 0;
			idx[i] = i;
		}//初始化
		for (int j = 0; j < nobj; j++)
		{
			for (int i = 0; i < s2; i++)
			{
				idd[i] = i;
				obj[i] = list[i].y_obj[j];
			}
			minfastsort(obj, idd, s2, s2);
			density[idd[0]] += -1.0e+30;//边界值取无穷大
			density[idd[s2 - 1]] += -1.0e+30;
			for (int k = 1; k < s2 - 1; k++)
			{
				density[idd[k]] += -(obj[k] - obj[k - 1] + obj[k + 1] - obj[k]);//求各个点的拥挤度距离								
				//cout<<"density[] = "<<density[idd[k]]<<"\n\n";
			}
		}
		delete[] idd;
		delete[] obj;
		int s3 = pops - population.size();
		//	cout<<"s3 = "<<s3<<"\n\n";
		minfastsort(density, idx, s2, s3);
		for (int i = 0; i < s3; i++)
			population.push_back(list[idx[i]]);
		delete[] density;
		delete[] idx;
	}
	offspring.clear();
}

//将N个个体放到offspring中，并保证没有重复
void zc_NSGA2::fill_union(zc_NSGA2Ind &ind)
{
	bool flag = true;
	int size = offspring.size();//计算offspring中个体的数量
	//cout<<"计算offspring中个体的数量: "<<size<<"\n\n";
	for (int i = 0; i < size; i++)
		if (ind == offspring[i])
		{
			flag = false;
			break;
		}
	if (flag)
		offspring.push_back(ind);
}

//将数组按照从小到大排序
void minfastsort(double x[], int idx[], int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = i + 1; j < n; j++)
			if (x[i] > x[j])
			{
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				int temp2 = idx[i];
				idx[i] = idx[j];
				idx[j] = temp2;
			}
	}
}

/*static unsigned int my_int = 1;
double rnd_uni(long *l)
{
	//srand((unsigned int)*l);
	//srand(time(NULL));
	my_int = (my_int + 500) % (INT_MAX - 50);
	srand(my_int);
	return rand() * 1.0 / RAND_MAX;
}*/

union FloatRand
{
	struct
	{
		unsigned long Frac : 23;
		unsigned long Exp : 8;
		unsigned long Signed : 1;
	} BitArea;
	float Value;
	unsigned long Binary; /* for debug only */
};

double rnd_uni(long *l)
{
	/*union FloatRand r;

	r.BitArea.Signed = 0;
	r.BitArea.Exp = 1;
	r.BitArea.Frac = (rand() * rand()) % 0x800000;
	if (r.BitArea.Frac == 0x7FFFFF)
		r.BitArea.Exp = 0x7D;
	else if (r.BitArea.Frac == 0)
		r.BitArea.Exp = 0x7E;
	else
		r.BitArea.Exp = 0x7E;

	//return r.Value;
	return (double)(r.Value - 0.5)*2.0;*/

	return rand() * 1.0 / RAND_MAX;
}